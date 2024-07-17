using trajectory_optimization_game_theory 

using LinearAlgebra: norm 
using Statistics: mean 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister


## ============================================ ##
# init params 

# start time 
tt = 0 

# init game 
game = game_struct( [], [], [] ) 
push!( game.tt, tt ) 

# orbit parameters of pursuer and evader 
tof = 1000          # tof for pursuer to catch up to evader  
N   = 10            # segments 
mu  = 398600.4415   # gravitational parameter 
r   = 6378.0        # Earth radius 
params = ( mu = mu, r = r, N = N, tof = tof ) 

kep0_P = [ r+420.0, 0.1, 20*pi/180, 10.0*pi/180, 20.0*pi/180, 20.0*pi/180 ]
kep0_E = [ r+520.0, 0.1, 20*pi/180, 10.0*pi/180, 20.0*pi/180, 25.0*pi/180 ]
rv_0_E = kep2cart(kep0_E, mu) 
rv_0_P = kep2cart(kep0_P, mu) 

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_E = vv2m(rv_E) 
rv_P = vv2m(rv_P) 

# save reference orbit 
ref_orbit = copy( kep0_E ) 

# compute vertices of polygon 
rv_vec = rv_E[end,:] 
vertices = polygon_vertices( rv_vec ) 


## ============================================ ## 
# plotting stuff 

fig = plot_axes3d( ) 
fig = plot_orbit( rv_E, fig ) 
fig = plot_orbit( rv_P, fig ) 

# plot 
fig = plot_axes3d( )
fig = plot_orbit( rv_E, fig ) 
fig = plot_orbit( rv_P, fig ) 
fig = plot_polygon( rv_vec, fig ) 


## ============================================ ## 
# compute all possible Δv solutions 

# function all_Δv_solns(  ) 

# save player state and control hists 
p = player_struct( [], [], [], [], [], [], [] ) 
players = [ p, deepcopy(p) ]  

players[1].rv_0 = rv_E 
players[2].rv_0 = rv_P 

# PLAYER ONE (EVADER) 
players = player_XU( params, vertices, players ) 

# compute cost matrices 
players = player_cost_matrices( players ) 

# solve mixed nash 
players = choose_weights( players, rng )


## ============================================ ##

push!( game.p1_state, players[1] ) 
push!( game.p2_state, players[2] )  

## ============================================ ##

# now that chosen has been sampled, choose the players[ii].X_vertices, players[ii].U_vertices, and players[ii].t_vertices 

# propagate time and k_idx forward 
tt_replan = 5       # replan every 5 * tof/N (100) seconds!!! 
tt = game.tt[end] + tt_replan * params.tof / params.N 
push!( game.tt, tt ) 
k_replan = 1 


# propagate SC state forward 
p1 = game.p1_state[ k_replan ] 
p2 = game.p2_state[ k_replan ] 

rv_0_E = p1.X[ p1.chosen ][ tt_replan + 1, : ]
rv_0_P = p2.X[ p2.chosen ][ tt_replan + 1, : ]

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

# compute vertices of polygon 
rv_vec   = rv_E[end,:] 
vertices = polygon_vertices( rv_vec ) 

## ============================================ ##





