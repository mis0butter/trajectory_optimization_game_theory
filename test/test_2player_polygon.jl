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
params, players, game = init_game( rng, "mixed" ) ; 

# ----------------------- #
# next game steps  

N_replan = 10   
for ii = 1 : N_replan - 1 
    print( "step: ", ii, "\n" ) 
    game = prop_game_step( game, params, rng ) 
end 

# ----------------------- #
# plotting stuff 

# k = 1 
fig = plot_p1_p2_traj( game, params, N_replan )  


## ============================================ ##
## ============================================ ##
# run multiple games 

N_games  = 100 
N_replan = 10 

# games_vec = run_MC_games( rng, N_games, N_replan, "mixed" ) 
games_vec = run_MC_games( rng, N_games, N_replan, "pure" ) 
games_vec = run_MC_games( rng, N_games, N_replan, "random" ) 

















## ============================================ ##
# save 

save_folder   = string( "test/results/", params.strategy, "/" ) 
filename      = string( "games_", N_games, ".jld2" ) 
full_filename = string( save_folder, filename ) 

@save full_filename games_vec params 



## ============================================ ##
## ============================================ ##
## ============================================ ##


ii   = 9 
game = games_vec[ ii ] 
fig  = plot_game_stats( game, params ) 


## ============================================ ##

fig = plot_MC_stats( games_vec ) 


## ============================================ ##
# save as CSV 



## ============================================ ##
# read CSV 

ii = 1 
game_folder = string(save_folder, "game_", ii,"/") 
filename    = string(game_folder, "tt_hist.csv")

tt_hist     = CSV.read( filename, DataFrame ) 




