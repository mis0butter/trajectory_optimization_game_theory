using trajectory_optimization_game_theory 

using LinearAlgebra: norm 
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

# propagate SC state forward 
p1 = game.p1_state[ end ] 
p2 = game.p2_state[ end ] 

# get current state 
rv_E = p1.X[ p1.chosen ][ params.tt_replan + 1, : ]
rv_P = p2.X[ p2.chosen ][ params.tt_replan + 1, : ]

t_E, rv_E_hist = propagate_2Body(rv_E, tof, mu, 1.0) 
t_P, rv_P_hist = propagate_2Body(rv_P, tof, mu, 1.0) 
rv_P_hist = vv2m(rv_P_hist) 
rv_E_hist = vv2m(rv_E_hist) 

# compute rv_ref for vertices of polygon - NEEDS TO BE UPDATED!!!! 
rv_ref_polygon = rv_E_hist[end,:] 

# save player state and control hists 
p = player_struct( [], [], [], [], [], [], [] ) 
players = [ p, deepcopy(p) ]  

players[1].rv_0_hist = rv_E_hist 
players[2].rv_0_hist = rv_P_hist 

# propagate time and k_idx forward 
push!( game.tt, game.tt[end] + params.tt_step ) 
push!( game.k_replan, game.k_replan[end] + 1 )  
push!( game.rv_E, rv_E ) 
push!( game.rv_P, rv_P ) 
push!( game.rv_ref_E, rv_E ) # THIS NEEDS TO BE UPDATED TO REFERENCE ORBIT !!! 
push!( game.rv_ref_polygon, rv_ref_polygon )  

## ============================================ ##

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

rv_E_hist = game.p1_state[k].rv_0_hist 
rv_P_hist = game.p2_state[k].rv_0_hist 

# propagate SC state forward 
p1 = game.p1_state[ k ] 
p2 = game.p2_state[ k ] 

# get current state 
rv_E = p1.X[ p1.chosen ][ params.tt_replan + 1, : ]
rv_P = p2.X[ p2.chosen ][ params.tt_replan + 1, : ]

# plot 
fig = plot_axes3d(  ) 
# fig = plot_orbit( rv_E_hist, fig ) 
# fig = plot_orbit( rv_P_hist, fig ) 
fig = plot_polygon( game.rv_ref_polygon[k], fig ) 
fig = plot_Δv_weights( game, params, k, fig )      

for i in 1 : params.tt_replan 

    # plot player 1 
    rv_E = p1.X[ p1.chosen ][ i, : ] 
    fig2 = plot_scatter3d( rv_E[1], rv_E[2], rv_E[3], fig, :circle, :blue, 20 ) 

    # plot player 2 
    rv_P = p2.X[ p2.chosen ][ i, : ] 
    fig2 = plot_scatter3d( rv_P[1], rv_P[2], rv_P[3], fig2, :circle, :red, 20 ) 

end 

## ============================================ ##

k = 1 
fig = plot_p1_p2_traj( game, params, k )  

