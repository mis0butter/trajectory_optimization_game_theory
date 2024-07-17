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

# now that chosen has been sampled, choose the players[ii].X_vertices, players[ii].U_vertices, and players[ii].t_vertices 

# propagate time and k_idx forward 
tt = game.tt[end] + params.tt_step 
push!( game.tt, tt ) 
push!( game.k_replan, game.k_replan[end] + 1 )  

# propagate SC state forward 
p1 = game.p1_state[ game.k_replan[ end - 1 ] ] 
p2 = game.p2_state[ game.k_replan[ end - 1 ] ] 

rv_0_E = p1.X[ p1.chosen ][ params.tt_replan + 1, : ]
rv_0_P = p2.X[ p2.chosen ][ params.tt_replan + 1, : ]

t_E, rv_E_hist = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P_hist = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P_hist = vv2m(rv_P_hist) 
rv_E_hist = vv2m(rv_E_hist) 

# compute vertices of polygon 
rv_vec   = rv_E_hist[end,:] 
vertices = polygon_vertices( rv_vec ) 

## ============================================ ##




## ============================================ ## 
# plotting stuff 

# fig = plot_axes3d(  ) 
# fig = plot_orbit( rv_E_hist, fig ) 
# fig = plot_orbit( rv_P_hist, fig ) 

# # plot 
# fig = plot_axes3d(  )
# fig = plot_orbit( rv_E_hist, fig ) 
# fig = plot_orbit( rv_P_hist, fig ) 
# fig = plot_polygon( rv_ref_polygon, fig ) 


