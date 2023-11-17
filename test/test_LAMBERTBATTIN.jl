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
# plot 

using GLMakie 

text_offset = (0,10) 

# initialize figure 
fig = Figure() 
Axis3(fig[1, 1], 
    xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)", 
    title = "Lambert Solution") 
    
# plot lambert solution 
lines!( x_lambert[:,1], x_lambert[:,2], x_lambert[:,3]; linewidth = 2, label = "lambert" ) 
scatter!( x_lambert[1,1], x_lambert[1,2], x_lambert[1,3]; marker = :circle, markersize = 10, color = :black ) 
text!( x_lambert[1,1], x_lambert[1,2], x_lambert[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 

# plot reverse 
lines!( x_reverse[:,1], x_reverse[:,2], x_reverse[:,3]; linewidth = 2, label = "reverse" ) 
scatter!( x_reverse[1,1], x_reverse[1,2], x_reverse[1,3]; marker = :circle, markersize = 10, color = :black ) 
# text!( x_reverse[1,1], x_reverse[1,2], x_reverse[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 

# plot E 
lines!( x_E[:,1], x_E[:,2], x_E[:,3]; linewidth = 2, label = "reverse" ) 
scatter!( x_E[1,1], x_E[1,2], x_E[1,3]; marker = :circle, markersize = 10, color = :black ) 
text!( x_E[1,1], x_E[1,2], x_E[1,3]; text = "E IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 

# plot target 
scatter!( r2[1], r2[2], r2[3]; marker = :circle, markersize = 10, color = :black ) 
text!( r2[1], r2[2], r2[3]; text = "target", color = :gray, offset = text_offset, align = (:center, :bottom) )

display(fig) 



