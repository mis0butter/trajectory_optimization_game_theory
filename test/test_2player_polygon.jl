using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean, var, std 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

using CSV, DataFrames 
using Infiltrator 

rng = MersenneTwister( 1 ) 


## ============================================ ##
# init params 

rng = MersenneTwister( 1 ) 

# run_game( rng, N_replan = 10, p1_strategy = "mixed", p2_strategy = "mixed" ) 

N_replan = 2  
p1_strategy  = "mixed" 
p2_strategy  = "mixed" 
game, params = run_game( rng, N_replan, p1_strategy, p2_strategy ) 

# ----------------------- #
# plotting stuff 

# k = 1 
# fig = plot_p1_p2_traj( game, params, N_replan )  
fig = plot_game_stats( game, params ) 



## ============================================ ## 
## ============================================ ## 
# run multiple games 

N_games  = 10 
N_replan = 10 

games_vec = run_MC_games( rng, N_games, N_replan, "mixed" ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "pure" ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "random" ) 

fig = plot_MC_stats( games_vec ) 


## ============================================ ##
# load and plot games_vec 

games_vec = load_games_vec( 100, "mixed" ) 
fig = plot_MC_stats( games_vec ) 


















## ============================================ ##
## ============================================ ##
# plot single game 

ii   = 9 
game = games_vec[ ii ] 
fig  = plot_game_stats( game, params ) 




## ============================================ ##
# read CSV 

ii = 1 
game_folder = string(save_folder, "game_", ii,"/") 
filename    = string(game_folder, "tt_hist.csv")

tt_hist     = CSV.read( filename, DataFrame ) 




