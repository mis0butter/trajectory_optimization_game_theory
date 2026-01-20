using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean, var, std 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

using CSV, DataFrames 
using Infiltrator 


## ====================================================================
## run single game 
## ====================================================================

rng = MersenneTwister( 1 ) 

# run_game( rng, k_replan = 10, p1_strategy = "mixed", p2_strategy = "mixed" ) 

k_replan = 2 
p1_strategy  = 1 
p2_strategy  = 1
game, params = run_game( rng, k_replan, p1_strategy, p2_strategy ) 

# ---------------------------------- 
# plotting stuff 

# k = 1 
fig = plot_p1_p2_traj( game, params, k_replan )  
fig_stats = plot_game_stats( game, params ) 


## ==================================================================== 
## test running multiple games 
## ==================================================================== 

rng = MersenneTwister( 1 ) 

N_games  = 100 
k_replan = 10 

p1_strategy  = "greedy" 
p2_strategy  = "mixed" 

# games_vec = run_MC_games( rng, N_games, k_replan, p1_strategy, p2_strategy ) 
# games_vec = run_MC_games(   rng, N_games, k_replan, "random",    "greedy"    ) 
# games_vec = run_MC_games(   rng, N_games, k_replan, "random",    "mixed"   ) 
# games_vec = run_MC_games(   rng, N_games, k_replan, "mixed",     "greedy"    ) 
# games_vec = run_MC_games( rng, N_games, k_replan, "greedy" ) 
# games_vec = run_MC_games( rng, N_games, k_replan, "random" ) 

games_vec = run_MC_games_parallel( rng, N_games, k_replan, p1_strategy, p2_strategy ) 

## ==================================================================== 

fig = plot_MC_stats( games_vec ) 


