using trajectory_optimization_game_theory

using LinearAlgebra: norm
using Statistics: mean, std
using Printf: @sprintf

using Random: MersenneTwister
using CSV, DataFrames 

## ==================================================================== 
## Belief convergence helpers 
## ==================================================================== 

function belief_entropy(belief::Vector{Float64})
    prob = belief ./ sum(belief)
    -sum(p_i * log(max(p_i, 1e-30)) for p_i in prob)
end

# ----------------------------- 

function correct_id_rate(games_vec, player_idx, opponent_strategy)
    tracked = games_vec[1].p1_state[end].tracked_strategies
    true_idx = findfirst(s -> s == opponent_strategy, tracked)
    if isnothing(true_idx)
        return NaN
    end
    n_correct = 0
    for game in games_vec
        state = player_idx == 1 ? game.p1_state[end] : game.p2_state[end]
        if argmax(state.strategy_belief) == true_idx
            n_correct += 1
        end
    end
    return n_correct / length(games_vec)
end

# ----------------------------- 

function mean_belief_entropy(games_vec, player_idx)
    entropies = Float64[]
    for game in games_vec
        state = player_idx == 1 ? game.p1_state[end] : game.p2_state[end]
        push!(entropies, belief_entropy(Float64.(state.strategy_belief)))
    end
    return mean(entropies)
end 

## ==================================================================== 
## Run all MC cases across multiple orbital periods 
## ==================================================================== 

rng = MersenneTwister(1)

N_games = 100
n_game_steps_vec = [10, 20, 30, 50]

strategies = ["greedy", "mixed", "random", "FP_greedy", "FP_mixed", "Meta_greedy", "Meta_mixed"]

results_df = DataFrame(
    n_game_steps=Int[],
    p1_strategy=String[],
    p2_strategy=String[],
    game_value=Float64[],
    game_value_last_quarter=Float64[],
    p1_U_total=Float64[],
    p2_U_total=Float64[],
    player_distance=Float64[],
    player_distance_last_quarter=Float64[],
    p1_belief_entropy=Float64[],
    p2_belief_entropy=Float64[],
    p1_correct_id_rate=Float64[],
    p2_correct_id_rate=Float64[],
)

# ----------------------------- 

for n_game_steps in n_game_steps_vec

    n_orbits = n_game_steps / 10
    println("\n" * "="^80)
    println("n_game_steps = $n_game_steps  ($n_orbits orbits)")
    println("="^80)

    for p1_strategy in strategies
        for p2_strategy in strategies

            println("Running: p1 = $p1_strategy, p2 = $p2_strategy ...")

            games_vec = run_MC_games_parallel(rng, N_games, n_game_steps, p1_strategy, p2_strategy)

            # ---- compute stats ---- 
            local_params = games_vec[1].params[1]
            sprintf_stats = print_MC_stats(games_vec, local_params, false)

            stage_cost_games, games_costs_mean, _ = stage_cost_games_fn(games_vec)
            game_val = mean(games_costs_mean)

            # ---- late-game metrics (last 25% of time steps) ---- 
            n_t = length(games_costs_mean)
            i_lq = max(1, n_t - n_t ÷ 4 + 1)
            game_val_lq = mean(stage_cost_games[:, i_lq:end])

            stats = MC_stats(games_vec, local_params)
            n_d = length(stats.dist_norm_mean)
            i_lq_d = max(1, n_d - n_d ÷ 4 + 1)
            dist_lq = mean(stats.dist_norm_mean[i_lq_d:end])

            # ---- belief convergence diagnostics ---- 
            p1_ent = mean_belief_entropy(games_vec, 1)
            p2_ent = mean_belief_entropy(games_vec, 2)
            p1_cid = correct_id_rate(games_vec, 1, p2_strategy)
            p2_cid = correct_id_rate(games_vec, 2, p1_strategy)

            push!(results_df, (
                n_game_steps=n_game_steps,
                p1_strategy=p1_strategy,
                p2_strategy=p2_strategy,
                game_value=game_val,
                game_value_last_quarter=game_val_lq,
                p1_U_total=parse(Float64, sprintf_stats.p1_Unorm_mean_end),
                p2_U_total=parse(Float64, sprintf_stats.p2_Unorm_mean_end),
                player_distance=parse(Float64, sprintf_stats.dist_norm_mean_mean),
                player_distance_last_quarter=dist_lq,
                p1_belief_entropy=p1_ent,
                p2_belief_entropy=p2_ent,
                p1_correct_id_rate=p1_cid,
                p2_correct_id_rate=p2_cid,
            ))

            println("  done. game_value = $(@sprintf("%.3g", game_val)), " *
                    "game_val_lq = $(@sprintf("%.3g", game_val_lq)), " *
                    "dist = $(sprintf_stats.dist_norm_mean_mean), " *
                    "dist_lq = $(@sprintf("%.3g", dist_lq)), " *
                    "p1_ent = $(@sprintf("%.3f", p1_ent)), " *
                    "p2_ent = $(@sprintf("%.3f", p2_ent))")
        end
    end
end

## ==================================================================== 
## Save combined multi-orbit results as CSV 
## ==================================================================== 

CSV.write("test/MC_results_multi_orbit.csv", results_df)
println("\nSaved multi-orbit results to test/MC_results_multi_orbit.csv") 

# ## ==================================================================== 
# ## Print pivot tables (like the image) 
# ## ==================================================================== 

# println("\n" * "="^80)
# println("GAME VALUE TABLE")
# println("="^80)

# # header 
# print(rpad("", 15))
# for p1 in strategies
#     print(rpad("p1 $p1", 15))
# end
# println()

# for p2 in strategies
#     print(rpad("p2 $p2", 15))
#     for p1 in strategies
#         row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
#         val = @sprintf "%.3g" row.game_value[1]
#         print(rpad(val, 15))
#     end
#     println()
# end


# println("\n" * "="^80)
# println("P1 U TOTAL TABLE")
# println("="^80)

# print(rpad("", 15))
# for p1 in strategies
#     print(rpad("p1 $p1", 15))
# end
# println()

# for p2 in strategies
#     print(rpad("p2 $p2", 15))
#     for p1 in strategies
#         row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
#         val = @sprintf "%.3g" row.p1_U_total[1]
#         print(rpad(val, 15))
#     end
#     println()
# end


# println("\n" * "="^80)
# println("P2 U TOTAL TABLE")
# println("="^80)

# print(rpad("", 15))
# for p1 in strategies
#     print(rpad("p1 $p1", 15))
# end
# println()

# for p2 in strategies
#     print(rpad("p2 $p2", 15))
#     for p1 in strategies
#         row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
#         val = @sprintf "%.3g" row.p2_U_total[1]
#         print(rpad(val, 15))
#     end
#     println()
# end


# println("\n" * "="^80)
# println("PLAYER DISTANCE TABLE")
# println("="^80)

# print(rpad("", 15))
# for p1 in strategies
#     print(rpad("p1 $p1", 15))
# end
# println()

# for p2 in strategies
#     print(rpad("p2 $p2", 15))
#     for p1 in strategies
#         row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
#         val = @sprintf "%.3g" row.player_distance[1]
#         print(rpad(val, 15))
#     end
#     println()
# end

# println("\nDone!")
