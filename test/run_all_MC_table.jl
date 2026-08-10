using trajectory_optimization_game_theory

using Random: MersenneTwister

## ==================================================================== 
## Run all MC cases — only n=30 needed since analyze_late_game.jl
## truncates to windows [10, 20, 30] from a single run.
## ==================================================================== 

rng = MersenneTwister(1)

N_games = 50
n_game_steps = 30

strategies = ["greedy", "mixed", "random", "FP_greedy", "FP_mixed"]
# Full set if needed later:
# strategies = ["greedy", "mixed", "random", "FP_greedy", "FP_mixed", "Meta_greedy", "Meta_mixed"]

n_total = length(strategies)^2

global i_run = 0

for p1_strategy in strategies
    for p2_strategy in strategies

        global i_run += 1
        println("\n[$i_run/$n_total] p1=$p1_strategy vs p2=$p2_strategy  (n_game_steps=$n_game_steps)")

        games_vec = run_MC_games_parallel(rng, N_games, n_game_steps, p1_strategy, p2_strategy)

        println("  saved $(length(games_vec)) games")
    end
end

println("\nDone — run analyze_late_game.jl to produce CSV.")
