using trajectory_optimization_game_theory 

using LinearAlgebra: norm 
using Statistics: mean 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister


## ============================================ ##
# init params 

# set up rng 
rng = MersenneTwister(1) 

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
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

## ============================================ ## 
# plotting stuff 

# compute vertices of polygon 
rv_vec = rv_E[end,:] 

vertices = polygon_vertices( rv_vec ) 

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
players = [] 

# PLAYER ONE (EVADER) 
players, fig = player_XU( params, rv_E, vertices, players, fig ) 

# PLAYER TWO (PURSUER) 
players, fig = player_XU( params, rv_P, vertices, players, fig ) 

# compute cost matrices 
player1_cost_matrix, player2_cost_matrix = player_cost_matrices( players ) 
players[1].cost = player1_cost_matrix   
players[2].cost = player2_cost_matrix   

## ============================================ ##
# solve mixed nash 

# mixing weights - ZERO SUM GAME!!! 
mixing_weights = let
    sol = solve_mixed_nash( players[1].cost )
    (; sol.x, sol.y) 
end 
players[1].weights = mixing_weights[1] 
players[2].weights = mixing_weights[2] 

# sample from mixed nash 
chosen = [sample(rng, ProbabilityWeights(weights)) for weights in mixing_weights] 
players[1].chosen = chosen[1] 
players[2].chosen = chosen[2] 


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
p1 = game.p1_state[k_replan] 
p2 = game.p2_state[k_replan] 

rv_0_E = p1.X[ p1.chosen ][ tt_replan + 1, : ]
rv_0_P = p2.X[ p2.chosen ][ tt_replan + 1, : ]

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

# compute vertices of polygon 
rv_vec   = rv_E[end,:] 
vertices = polygon_vertices( rv_vec ) 

# save player state and control hists 
players = [] 

# PLAYER ONE (EVADER) 
players, fig = player_XU( params, rv_E, vertices, players, fig ) 

# PLAYER TWO (PURSUER) 
players, fig = player_XU( params, rv_P, vertices, players, fig ) 

# compute cost matrices 
player1_cost_matrix, player2_cost_matrix = player_cost_matrices( players ) 
cost_matrices = ( player1_cost_matrix, player2_cost_matrix ) 

# mixing weights 
mixing_weights = let
    sol = solve_mixed_nash( cost_matrices[1] )
    (; sol.x, sol.y)
end 

# sample from mixed nash 
chosen = [sample(rng, ProbabilityWeights(weights)) for weights in mixing_weights]

## ============================================ ## 



