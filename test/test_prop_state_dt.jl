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

rv0 = [ r1; v1 ] 
t_lambert, rv_hist = propagate_2Body( rv0, tof, mu, 1.0 ) 
rv_hist = mapreduce( permutedims, vcat, rv_hist ) 

## ============================================ ##
# break up delta v into smaller segments 

# initial position has all z velocity 
vz  = [ 0; 0; norm(v1) ]

# compute delta v vec for lambert solution 
dv  = v1 - vz   

# break up delta v into smaller segments 
N = 10 
dv_vec = [] 
for i = 1 : N 
    push!( dv_vec, dv/N ) 
end 
dv_vec = mapreduce( permutedims, vcat, dv_vec ) 

## ============================================ ##
# propagate each segment 

dt   = tof / N 
i    = 1 
xk   = rv0 

rv_hist = [] 
for i = 1 : N 
    xkdv = apply_dv( xk, dv_vec[i,:]) 
    t, x = propagate_2Body( xkdv, dt, mu ) 
    for j = 1 : length(x) 
        push!( rv_hist, x[j] ) 
    end 
    xk = x[end]
end 
rv_hist = mapreduce( permutedims, vcat, rv_hist ) 

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
lines!( rv_hist[:,1], rv_hist[:,2], rv_hist[:,3]; linewidth = 2, label = "lambert" ) 
scatter!( rv_hist[1,1], rv_hist[1,2], rv_hist[1,3]; marker = :circle, markersize = 10, color = :black ) 
text!( rv_hist[1,1], rv_hist[1,2], rv_hist[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 

display(fig) 
