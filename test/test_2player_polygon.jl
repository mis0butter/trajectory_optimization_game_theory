using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

rng = MersenneTwister( 1 )


## ============================================ ##
# init params 

params, players, game = init_game(  ) 

# compute all possible Δv solutions 
players = players_states( params, game, players, rng ) 

# save player state in game 
push!( game.p1_state, players[1] ) 
push!( game.p2_state, players[2] )  


## ============================================ ##
# set up for next game step 

game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 


## ============================================ ##
## ============================================ ## 
# plotting stuff 

# k_replan step of the game 
k = 5 
fig = plot_p1_p2_traj( game, params, k )  

