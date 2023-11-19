using trajectory_optimization_game_theory
using LinearAlgebra 

## ============================================ ## 
# define IC, target state, and lambert solve 

r1      = [20.0e6, 20.0e6, 0]   # [m] 
r2      = [-20.0e6, 10.0e6, 0]  # [m] 
tof     = 1.0 * 86400 
mu      = 398600.4418e9         # [m^3/s^2] 
dm      = "pro" 
Dtsec   = tof 

## ============================================ ##
# solve and propagate lambert orbit 

v1, v2 = lambertbattin( r1, r2, mu, dm, tof ) 

x0 = [ r1; v1 ] 
t_lambert, x_lambert = propagate_2Body( x0, tof, mu, 1.0 ) 
x_lambert = mapreduce( permutedims, vcat, x_lambert ) 

## ============================================ ##
# test negative propagation 

# xf = x_lambert[end,:] 
vf = v2 ; vf[end] = 100
# xf = [r2; v2]
xf = [r2; vf] 
t_reverse, x_reverse = propagate_2Body( xf, -tof, mu, 1.0 ) 
x_reverse = mapreduce( permutedims, vcat, x_reverse ) 

x_E_0  = x_reverse[end,:] 
t_E, x_E = propagate_2Body( x_E_0, tof, mu, 1.0 ) 
x_E = mapreduce( permutedims, vcat, x_E ) 

## ============================================ ##

fig = plot_orbit( rv_lambert )
fig = plot_orbit( x_reverse, fig )  
fig = plot_orbit( x_E, fig )  

# plot target 
scatter!( r2[1], r2[2], r2[3]; marker = :circle, markersize = 10, color = :black ) 
text!( r2[1], r2[2], r2[3]; text = "target", color = :gray, offset = (0,-10), align = (:center, :bottom) )




