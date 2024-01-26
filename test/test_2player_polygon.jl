using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Statistics 
using Optim 

## ============================================ ##
# init params 

mu = 398600.4415
r  = 6378.0
kep0_P = [ r+420.0, 0.1, 20*pi/180, 10.0*pi/180, 20.0*pi/180, 20.0*pi/180 ]
rv_0_P = kep2cart(kep0_P, mu) 
kep0_E = [ r+450.0, 0.2, 10.6*pi/180, 12.0*pi/180, 0.0, 50.0*pi/180 ]
rv_0_E = kep2cart(kep0_E, mu) 

# tof for pursuer to catch up to evader 
tof = 1000 

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

rv_vec = rv_E[end,:] 

# get vertices of polygon 
vertices = polygon_vertices( rv_vec ) 

fig = plot_axes3d( )
fig = plot_orbit( rv_E, fig ) 
fig = plot_orbit( rv_P, fig ) 

## ============================================ ##

# plot 
fig = plot_axes3d( )
fig = plot_orbit( rv_E, fig ) 
# fig = plot_orbit( rv_P, fig ) 
fig = plot_polygon( rv_vec, fig ) 

# segments 
N = 10 

# save player state and control hists 
players = [] 

# ----------------------- #
# PLAYER ONE 

# init state and end velocity (probably doesn't matter) 
rv_0 = rv_0_E 
v_f  = rv_E[end,4:6] 
# Δv_sol = min_Δv( rv_0, rv_f, tof, N, mu ) 

X_vertices = [] 
U_vertices = [] 
t_vertices = [] 
for i in eachindex(vertices)  

    rv_f   = [ vertices[i] ; v_f ]  
    Δv_sol = min_Δv_dist( rv_0, rv_f, tof, N, mu ) 
    t, rv_hist = prop_kepler_tof_Nseg( rv_0, Δv_sol, N, tof / N, mu ) 

    # save hist 
    push!( X_vertices, rv_hist ) 
    push!( U_vertices, Δv_sol ) 
    push!( t_vertices, t ) 

    fig    = plot_prop_Δv( rv_0, Δv_sol, N, tof / N, mu, fig )     

end 

push!( players, ( X_vertices = X_vertices, U_vertices = U_vertices, t_vertices = t_vertices ) ) 

# ----------------------- #
# PLAYER TWO 

# init state and end velocity (probably doesn't matter) 
rv_0 = rv_0_P  
v_f  = rv_E[end,4:6] 
# Δv_sol = min_Δv( rv_0, rv_f, tof, N, mu ) 

X_vertices = [] 
U_vertices = [] 
t_vertices = [] 
for i in eachindex(vertices)  

    rv_f   = [ vertices[i] ; v_f ]  
    Δv_sol = min_Δv_dist( rv_0, rv_f, tof, N, mu ) 
    t, rv_hist = prop_kepler_tof_Nseg( rv_0, Δv_sol, N, tof / N, mu ) 

    # save hist 
    push!( X_vertices, rv_hist ) 
    push!( U_vertices, Δv_sol ) 
    push!( t_vertices, t ) 

    fig    = plot_prop_Δv( rv_0, Δv_sol, N, tof / N, mu, fig )     

end 

push!( players, ( X_vertices = X_vertices, U_vertices = U_vertices, t_vertices = t_vertices ) ) 

## ============================================ ## 
# zero-sum game 

# game cost 
function stage_cost(x1, x2, u1, u2)
    sqrt(norm(x1[1:3] - x2[1:3]) + 0.1) + 0.1 * (norm(u1) - norm(u2))
end

# loop through time corresponding with control inputs 
us = eachindex( players[begin].U_vertices[begin][:,begin] )

# start with vertex 1 for players 1 and 2 
i_vert = 1 
j_vert = 1 


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
    cost2 = -stage_cost( x1, x2, u1, u2 )
    push!( player1_cost_tt, cost1 ) 
    push!( player2_cost_tt, cost2 ) 
    
end 
player1_cost = mean( player1_cost_tt )
player2_cost = mean( player2_cost_tt )  





# what do I need? 
# I need to save the state and control history for each player 
# I need to save the cost history for each player 

