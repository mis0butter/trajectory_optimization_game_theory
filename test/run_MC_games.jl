using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean, var, std 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

using CSV, DataFrames 
using Infiltrator 

rng = MersenneTwister( 1 ) 


## ====================================================================
## run single game 
## ====================================================================

rng = MersenneTwister( 1 ) 

# run_game( rng, N_replan = 10, p1_strategy = "mixed", p2_strategy = "mixed" ) 

N_replan = 10  
p1_strategy  = "mixed" 
p2_strategy  = "mixed" 
game, params = run_game( rng, N_replan, p1_strategy, p2_strategy ) 

# ---------------------------------- 
# plotting stuff 

# k = 1 
fig = plot_p1_p2_traj( game, params, N_replan )  
fig = plot_game_stats( game, params ) 


## ==================================================================== 
## test running multiple games 
## ==================================================================== 

rng = MersenneTwister( 1 ) 

N_games  = 100 
N_replan = 10 

p1_strategy  = "pure" 
p2_strategy  = "mixed" 

# games_vec = run_MC_games( rng, N_games, N_replan, p1_strategy, p2_strategy ) 
# games_vec = run_MC_games(   rng, N_games, N_replan, "random",    "pure"    ) 
# games_vec = run_MC_games(   rng, N_games, N_replan, "random",    "mixed"   ) 
# games_vec = run_MC_games(   rng, N_games, N_replan, "mixed",     "pure"    ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "pure" ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "random" ) 

games_vec = run_MC_games_parallel( rng, N_games, N_replan, p1_strategy, p2_strategy ) 

## ==================================================================== 

fig = plot_MC_stats( games_vec ) 


