using trajectory_optimization_game_theory

## ============================================ ## 
# define IC and target state 

mu = 398600.4415
r = 6378.0
kep0_P = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_E = [r+450.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]
t = (0.0, 3/4*orbitPeriod(kep0_E, mu)) 
prop_P = propagate_2Body(kep2cart(kep0_P, mu), t, mu)
prop_E = propagate_2Body(kep2cart(kep0_E, mu), t, mu)

x₀_P = prop_P.u[1] 
x₀_E = prop_E.u[1] 
xf_E = prop_E.u[end] 

## ============================================ ##
# lambert solve 

r1 = x₀_P[1:3] 
r2 = xf_E[1:3] 

tof     = t[end]
mu      = 398600.4418e9         # [m^3/s^2] 
dm      = "retro" 
Dtsec   = tof 
v1, v2  = lambertbattin(r1, r2, mu, dm, tof) 

x₀_P_lambert   = [r1; v1] 
prop_P_lambert = propagate_2Body( x₀_P_lambert, tof, mu )

## ============================================ ##
# massage data for plotting 

x_P = prop_P.u 
x_P = mapreduce( permutedims, vcat, x_P ) 

x_E = prop_E.u 
x_E = mapreduce( permutedims, vcat, x_E ) 

x_P_lambert = prop_P_lambert.u 
x_P_lambert = mapreduce( permutedims, vcat, x_P_lambert ) 

## ============================================ ##

using GLMakie 

# initialize figure 
fig = Figure() 
Axis3(fig[1, 1], 
    xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)", 
    title = "Lambert Solution") 

lines!( x_P[:,1], x_P[:,2], x_P[:,3]; linewidth = 2, label = "OG" ) 
lines!( x_E[:,1], x_E[:,2], x_E[:,3]; linewidth = 2, label = "Target" ) 
lines!( x_P_lambert[:,1], x_P_lambert[:,2], x_P_lambert[:,3]; linewidth = 2, label = "Lambert" ) 

display(fig) 











