using trajectory_optimization_game_theory

using JLD2
using LinearAlgebra: norm
using Statistics: mean
using Printf: @sprintf
using CSV, DataFrames


## ====================================================================
## Helper functions
## ====================================================================

function belief_entropy(belief::Vector{Float64})
    prob = belief ./ sum(belief)
    -sum(p_i * log(max(p_i, 1e-30)) for p_i in prob)
end

function compute_stage_costs(game, params, n_steps)
    n_seg = params.n_seg_per_game_step
    costs = Float64[]
    for kk = 1:n_steps
        p1 = game.p1_state[kk]
        p2 = game.p2_state[kk]
        for ii = 1:n_seg
            x1 = p1.X[p1.chosen][ii, :]
            u1 = p1.U[p1.chosen][ii, :]
            x2 = p2.X[p2.chosen][ii, :]
            u2 = p2.U[p2.chosen][ii, :]
            push!(costs, stage_cost(x1, x2, u1, u2))
        end
    end
    return costs
end

function compute_distances(game, params, n_steps)
    n_seg = params.n_seg_per_game_step
    dists = Float64[]
    for kk = 1:n_steps
        p1 = game.p1_state[kk]
        p2 = game.p2_state[kk]
        for ii = 1:n_seg
            r1 = p1.X[p1.chosen][ii, 1:3]
            r2 = p2.X[p2.chosen][ii, 1:3]
            push!(dists, norm(r1 - r2))
        end
    end
    return dists
end

function compute_U_total(game, params, n_steps, player_idx)
    n_seg = params.n_seg_per_game_step
    total = 0.0
    for kk = 1:n_steps
        p = player_idx == 1 ? game.p1_state[kk] : game.p2_state[kk]
        U = p.U[p.chosen][1:n_seg, :]
        for ii = 1:n_seg
            total += norm(U[ii, :])
        end
    end
    return total
end


## ====================================================================
## Post-process .jld2 files with truncated windows
## ====================================================================

results_dir = "test/results"
windows = [10, 20, 30]
# valid_strategies = Set(["greedy", "mixed", "random", "FP_greedy", "FP_mixed", "Meta_greedy", "Meta_mixed"])
valid_strategies = Set(["greedy", "mixed", "random", "FP_greedy", "FP_mixed"])

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

# Collect .jld2 paths from nested layout (test/results/n<k>/p1_X_p2_Y/)
jld2_paths = String[]
for entry in readdir(results_dir, join=true)
    isdir(entry) && startswith(basename(entry), "n") || continue
    for sub in readdir(entry, join=true)
        if isdir(sub) && startswith(basename(sub), "p1_")
            for f in readdir(sub, join=true)
                if endswith(f, ".jld2") && startswith(basename(f), "games_")
                    push!(jld2_paths, f)
                end
            end
        end
    end
end
sort!(jld2_paths)

for jld2_path in jld2_paths

    # Parse strategy names from parent folder: "p1_X_p2_Y"
    folder = basename(dirname(jld2_path))
    parts = split(folder, "_p2_")
    length(parts) == 2 || continue
    p1_strategy = String(parts[1][4:end])
    p2_strategy = String(parts[2])

    p1_strategy in valid_strategies && p2_strategy in valid_strategies || continue

    println("Loading $folder ...")
    @load jld2_path games_vec

    k_max = games_vec[1].i_game_step[end]
    params = games_vec[1].params[1]
    tracked = games_vec[1].p1_state[1].tracked_strategies

    p2_true_idx = findfirst(s -> s == p2_strategy, tracked)
    p1_true_idx = findfirst(s -> s == p1_strategy, tracked)

    for n_steps in windows

        if n_steps > k_max
            continue
        end

        game_costs_all = Float64[]
        game_costs_lq_all = Float64[]
        dist_all = Float64[]
        dist_lq_all = Float64[]
        p1_U_all = Float64[]
        p2_U_all = Float64[]
        p1_ent_all = Float64[]
        p2_ent_all = Float64[]
        p1_correct = 0
        p2_correct = 0

        for game in games_vec

            costs = compute_stage_costs(game, params, n_steps)
            push!(game_costs_all, mean(costs))
            n_c = length(costs)
            i_lq = max(1, n_c - n_c ÷ 4 + 1)
            push!(game_costs_lq_all, mean(costs[i_lq:end]))

            dists = compute_distances(game, params, n_steps)
            push!(dist_all, mean(dists))
            n_d = length(dists)
            i_lq_d = max(1, n_d - n_d ÷ 4 + 1)
            push!(dist_lq_all, mean(dists[i_lq_d:end]))

            push!(p1_U_all, compute_U_total(game, params, n_steps, 1))
            push!(p2_U_all, compute_U_total(game, params, n_steps, 2))

            p1_belief = Float64.(game.p1_state[n_steps].strategy_belief)
            p2_belief = Float64.(game.p2_state[n_steps].strategy_belief)
            push!(p1_ent_all, belief_entropy(p1_belief))
            push!(p2_ent_all, belief_entropy(p2_belief))

            if !isnothing(p2_true_idx) && argmax(p1_belief) == p2_true_idx
                p1_correct += 1
            end
            if !isnothing(p1_true_idx) && argmax(p2_belief) == p1_true_idx
                p2_correct += 1
            end
        end

        N = length(games_vec)

        push!(results_df, (
            n_game_steps=n_steps,
            p1_strategy=p1_strategy,
            p2_strategy=p2_strategy,
            game_value=mean(game_costs_all),
            game_value_last_quarter=mean(game_costs_lq_all),
            p1_U_total=mean(p1_U_all),
            p2_U_total=mean(p2_U_all),
            player_distance=mean(dist_all),
            player_distance_last_quarter=mean(dist_lq_all),
            p1_belief_entropy=mean(p1_ent_all),
            p2_belief_entropy=mean(p2_ent_all),
            p1_correct_id_rate=isnothing(p2_true_idx) ? NaN : p1_correct / N,
            p2_correct_id_rate=isnothing(p1_true_idx) ? NaN : p2_correct / N,
        ))
    end

    println("  Done: $p1_strategy vs $p2_strategy (k_max=$k_max)")
end


## ====================================================================
## Save results
## ====================================================================

sort!(results_df, [:n_game_steps, :p1_strategy, :p2_strategy])

CSV.write("test/MC_results_late_game.csv", results_df)
println("\nSaved results to test/MC_results_late_game.csv")
println("Total rows: $(nrow(results_df))")
