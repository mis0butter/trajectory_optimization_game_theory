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

# propagate current SC state forward 
p1 = game.p1_state[ end ] 
p2 = game.p2_state[ end ] 

# get current state 
rv_E = p1.X[ p1.chosen ][ params.tt_replan + 1, : ] 
rv_P = p2.X[ p2.chosen ][ params.tt_replan + 1, : ] 

rv_ref_E_polygon_hist = game.rv_ref_E_polygon[ end ] 

# find smallest angle between rv_E and rv_ref_E_polygon_hist and index 
cos_min = 100 
ii_min  = 1 
for ii in axes( rv_ref_E_polygon_hist, 1 )

    rv_ref_E_polygon = rv_ref_E_polygon_hist[ii,:] 
    dot_p = dot( rv_E, rv_ref_E_polygon ) / ( norm(rv_E) * norm(rv_ref_E_polygon) )  
    cos_a = acos( dot_p ) 

    if cos_a < cos_min  
        cos_min = cos_a  
        ii_min  = ii 
    end 

end 

# save reference orbit 
rv_ref_E  = rv_ref_E_polygon_hist[ii_min,:] 
kep_ref_E = cart2kep( rv_ref_E, params.mu ) 

# save OG reference orbit 
kep0_ref_E = params.kep0_ref_E 
kep0_ref_E[end] = kep_ref_E[end] 

rv0_ref_E = kep2cart( kep0_ref_E, params.mu ) 

t_E, rv_ref_E_polygon_hist = propagate_2Body(rv0_ref_E, params.tof, params.mu, 1.0) 
rv_ref_E_polygon_hist = vv2m(rv_ref_E_polygon_hist) 

# now move game forward one step 
push!( game.tt, game.tt[end] + params.tt_step ) 
push!( game.k_replan, game.k_replan[end] + 1 ) 
push!( game.rv_E, rv_E ) 
push!( game.rv_P, rv_P ) 
push!( game.rv_ref_E_polygon, rv_ref_E_polygon_hist ) 

# save player state and control hists 
p = player_struct( [], [], [], [], [], [], [] ) 
players = [ p, deepcopy(p) ]  

t_E, rv_E_hist = propagate_2Body(rv_E, params.tof, params.mu, 1.0) 
t_P, rv_P_hist = propagate_2Body(rv_P, params.tof, params.mu, 1.0) 
rv_P_hist = vv2m(rv_P_hist) 
rv_E_hist = vv2m(rv_E_hist) 

players[1].rv_0_hist = rv_E_hist 
players[2].rv_0_hist = rv_P_hist 

# compute all possible Δv solutions 
players = players_states( params, game, players, rng ) 

# save player state in game 
push!( game.p1_state, players[1] ) 
push!( game.p2_state, players[2] ) 


## ============================================ ##
## ============================================ ## 
# plotting stuff 

# k_replan step of the game 
k = 1 
fig = plot_p1_p2_traj( game, params, k )  

