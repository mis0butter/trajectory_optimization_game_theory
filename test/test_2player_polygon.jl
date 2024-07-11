using trajectory_optimization_game_theory 

using LinearAlgebra: norm 
using Statistics: mean 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister


## ============================================ ##
# init params 

rng = MersenneTwister(1) 

mu = 398600.4415
r  = 6378.0
kep0_P = [ r+420.0, 0.1, 20*pi/180, 10.0*pi/180, 20.0*pi/180, 20.0*pi/180 ]
rv_0_P = kep2cart(kep0_P, mu) 
kep0_E = [ r+520.0, 0.1, 20*pi/180, 10.0*pi/180, 20.0*pi/180, 25.0*pi/180 ]
rv_0_E = kep2cart(kep0_E, mu) 

# tof for pursuer to catch up to evader 
tof = 1000 

# segments 
N = 10 

params = ( mu = mu, r = r, N = N, tof = tof )

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

rv_vec = rv_E[end,:] 

# compute vertices of polygon 
vertices = polygon_vertices( rv_vec ) 

fig = plot_axes3d( )
fig = plot_orbit( rv_E, fig ) 
fig = plot_orbit( rv_P, fig ) 

## ============================================ ## 
# compute all possible Δv solutions 


# plot 
fig = plot_axes3d( )
fig = plot_orbit( rv_E, fig ) 
# fig = plot_orbit( rv_P, fig ) 
fig = plot_polygon( rv_vec, fig ) 

# function all_Δv_solns(  ) 

# save player state and control hists 
players = [] 

# PLAYER ONE (EVADER) 
players, fig = player_XU( params, rv_E, vertices, players, fig ) 

# PLAYER TWO (PURSUER) 
players, fig = player_XU( params, rv_P, vertices, players, fig ) 


## ============================================ ## 
# zero-sum game 

function player_cost_matrices( players ) 

    # game cost 
    function stage_cost(x1, x2, u1, u2)
        sqrt(norm(x1[1:3] - x2[1:3]) + 0.1) + 0.1 * (norm(u1) - norm(u2))
    end

    # loop through time corresponding with control inputs 
    us = eachindex( players[begin].U_vertices[begin][:,begin] )

    # start with vertex 1 for players 1 and 2 
    i_vert = 1 
    j_vert = 1 

    player1_cost_matrix = zeros(6, 6) 
    player2_cost_matrix = zeros(6, 6) 

    n_vertices = length( players[1].t_vertices ) 

    for i_vert in 1 : n_vertices 
        for j_vert in 1 : n_vertices 

            # loop through time 
            player1_cost_tt = [] 
            player2_cost_tt = [] 
            for tt in us 

                x1 = players[1].X_vertices[i_vert][tt,:] 
                u1 = players[1].U_vertices[i_vert][tt,:] 
                x2 = players[2].X_vertices[j_vert][tt,:] 
                u2 = players[2].U_vertices[j_vert][tt,:] 

                # compute costs for player 1 and 2  
                cost1 = stage_cost( x1, x2, u1, u2 )
                cost2 = - stage_cost( x1, x2, u1, u2 )
                push!( player1_cost_tt, cost1 ) 
                push!( player2_cost_tt, cost2 ) 
                
            end 
            player1_cost = mean( player1_cost_tt )
            player2_cost = mean( player2_cost_tt )  

            # save cost in matrix 
            player1_cost_matrix[i_vert, j_vert] = player1_cost 
            player2_cost_matrix[i_vert, j_vert] = player2_cost 

        end 
    end 

    return player1_cost_matrix, player2_cost_matrix 
end 

player1_cost_matrix, player2_cost_matrix = player_cost_matrices( players ) 



## ============================================ ##
# solve mixed nash 

cost_matrices = ( player1_cost_matrix, player2_cost_matrix ) 

# mixing weights 
mixing_weights = let
    sol = solve_mixed_nash( cost_matrices[1] )
    (; sol.x, sol.y)
end 
println( "mixing weights = ", mixing_weights ) 

# sample from mixed nash 
chosen = [sample(rng, ProbabilityWeights(weights)) for weights in mixing_weights]
println( "chosen = ", chosen ) 

## ============================================ ##

