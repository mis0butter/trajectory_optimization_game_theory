using trajectory_optimization_game_theory

using LinearAlgebra: norm
using Statistics: mean, std
using Printf: @sprintf

using Random: MersenneTwister
using CSV, DataFrames


## ==================================================================== 
## Run all MC cases and collect results into a table 
## ==================================================================== 

rng = MersenneTwister(1)

N_games = 100
k_replan = 10

strategies = ["greedy", "mixed", "random", "FP_greedy", "FP_mixed", "Meta_greedy", "Meta_mixed"]

# DataFrame to store results 
results_df = DataFrame(
    p1_strategy=String[],
    p2_strategy=String[],
    game_value=Float64[],
    p1_U_total=Float64[],
    p2_U_total=Float64[],
    player_distance=Float64[],
)

for p1_strategy in strategies
    for p2_strategy in strategies

        println("Running: p1 = $p1_strategy, p2 = $p2_strategy ...")

        games_vec = run_MC_games_parallel(rng, N_games, k_replan, p1_strategy, p2_strategy)

        # ---- compute stats ---- 
        sprintf_stats = print_MC_stats(games_vec, games_vec[1].params[1], false)

        # game value = mean of mean stage costs across all games 
        _, games_costs_mean, _ = stage_cost_games_fn(games_vec)
        game_val = mean(games_costs_mean)

        push!(results_df, (
            p1_strategy=p1_strategy,
            p2_strategy=p2_strategy,
            game_value=game_val,
            p1_U_total=parse(Float64, sprintf_stats.p1_Unorm_mean_end),
            p2_U_total=parse(Float64, sprintf_stats.p2_Unorm_mean_end),
            player_distance=parse(Float64, sprintf_stats.dist_norm_mean_mean),
        ))

        println("  done. game_value = $(@sprintf("%.3g", game_val)), " *
                "p1_U = $(sprintf_stats.p1_Unorm_mean_end), " *
                "p2_U = $(sprintf_stats.p2_Unorm_mean_end), " *
                "dist = $(sprintf_stats.dist_norm_mean_mean)")
    end
end


## ==================================================================== 
## Save full results as CSV 
## ==================================================================== 

CSV.write("test/MC_results_table.csv", results_df)
println("\nSaved full results to test/MC_results_table.csv")


## ==================================================================== 
## Print pivot tables (like the image) 
## ==================================================================== 

println("\n" * "="^80)
println("GAME VALUE TABLE")
println("="^80)

# header 
print(rpad("", 15))
for p1 in strategies
    print(rpad("p1 $p1", 15))
end
println()

for p2 in strategies
    print(rpad("p2 $p2", 15))
    for p1 in strategies
        row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
        val = @sprintf "%.3g" row.game_value[1]
        print(rpad(val, 15))
    end
    println()
end


println("\n" * "="^80)
println("P1 U TOTAL TABLE")
println("="^80)

print(rpad("", 15))
for p1 in strategies
    print(rpad("p1 $p1", 15))
end
println()

for p2 in strategies
    print(rpad("p2 $p2", 15))
    for p1 in strategies
        row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
        val = @sprintf "%.3g" row.p1_U_total[1]
        print(rpad(val, 15))
    end
    println()
end


println("\n" * "="^80)
println("P2 U TOTAL TABLE")
println("="^80)

print(rpad("", 15))
for p1 in strategies
    print(rpad("p1 $p1", 15))
end
println()

for p2 in strategies
    print(rpad("p2 $p2", 15))
    for p1 in strategies
        row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
        val = @sprintf "%.3g" row.p2_U_total[1]
        print(rpad(val, 15))
    end
    println()
end


println("\n" * "="^80)
println("PLAYER DISTANCE TABLE")
println("="^80)

print(rpad("", 15))
for p1 in strategies
    print(rpad("p1 $p1", 15))
end
println()

for p2 in strategies
    print(rpad("p2 $p2", 15))
    for p1 in strategies
        row = filter(r -> r.p1_strategy == p1 && r.p2_strategy == p2, results_df)
        val = @sprintf "%.3g" row.player_distance[1]
        print(rpad(val, 15))
    end
    println()
end

println("\nDone!")
