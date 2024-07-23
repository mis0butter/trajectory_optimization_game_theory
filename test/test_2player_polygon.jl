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

# # get most recent player states  
# rv_E, rv_P = rv_E_P( game, params ) 

# # return reference orbit for most recent step 
# rv_ref_E, kep_ref_E = find_ref_orbit( game, params ) 

# # generate rv_ref_E_polygon_hist 
# rv_ref_E_polygon_hist = ref_polygon_hist( kep_ref_E, params ) 

# # now move game forward one step 
# push!( game.tt, game.tt[end] + params.tt_step ) 
# push!( game.k_replan, game.k_replan[end] + 1 ) 
# push!( game.rv_E, rv_E ) 
# push!( game.rv_P, rv_P ) 
# push!( game.rv_ref_E_polygon, rv_ref_E_polygon_hist ) 

# # save player state and control hists 
# p = player_struct( [], [], [], [], [], [], [] ) 
# players = [ p, deepcopy(p) ]  

# _, rv_E_hist, _, rv_P_hist = prop_rv_E_P( rv_E, rv_P, params ) 

# players[1].rv_0_hist = rv_E_hist 
# players[2].rv_0_hist = rv_P_hist 

# # compute all possible Δv solutions 
# players = players_states( params, game, players, rng ) 

# # save player state in game 
# push!( game.p1_state, players[1] ) 
# push!( game.p2_state, players[2] ) 


## ============================================ ##
## ============================================ ## 
# plotting stuff 

# k_replan step of the game 
k = 1 
fig = plot_p1_p2_traj( game, params, k )  

