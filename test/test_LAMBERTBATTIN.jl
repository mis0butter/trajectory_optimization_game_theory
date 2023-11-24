using trajectory_optimization_game_theory

## ============================================ ## 
# define IC, target state, and lambert solve 

r1      = [20.0e6, 20.0e6, 0]   # [m] 
r2      = [-20.0e6, 10.0e6, 0]  # [m] 
tof     = 1.0 * 86400 
mu      = 398600.4418e9         # [m^3/s^2] 
dm      = "retro" 
Dtsec   = tof 

## ============================================ ##
# solve and propagate lambert orbit 

v1, v2  = lambertbattin(r1, r2, mu, dm, tof) 

x0 = [r1; v1] 
prop_P_lambert = propagate_2Body( x0, tof*1.2, mu, 1.0 ) 

x_lambert = prop_P_lambert.u 
x_lambert = mapreduce( permutedims, vcat, x_lambert ) 

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
lines!( x_lambert[:,1], x_lambert[:,2], x_lambert[:,3]; linewidth = 2, label = "OG" ) 
scatter!( x_lambert[1,1], x_lambert[1,2], x_lambert[1,3]; marker = :circle, markersize = 10, color = :black ) 
text!( x_lambert[1,1], x_lambert[1,2], x_lambert[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 

# plot target 
scatter!( r2[1], r2[2], r2[3]; marker = :circle, markersize = 10, color = :black ) 
text!( r2[1], r2[2], r2[3]; text = "target", color = :gray, offset = text_offset, align = (:center, :bottom) )

display(fig) 



