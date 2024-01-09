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
t_lambert, rv_lambert = propagate_2Body( rv0, tof, mu, 1.0 ) 
rv_lambert = mapreduce( permutedims, vcat, rv_lambert ) 

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
# prop_2Body_tof_Nseg 

# Creating Iteration Variables
xk = copy(rv0)
dt_N = tof / N 
Δt = N * dt_N 

# Propagating Through Each Segment 
X_hist = [ apply_Δv( xk, dv_vec[i,:] ) ] 
t_hist = [ 0 ] 
# for i = 1 : N 
i = 1 

    # apply dv 
    xkdv = apply_Δv( xk, dv_vec[i,:] ) 

    # propagate and save 
    t, x = propagate_2Body( xkdv, dt_N, mu ) 
    for j = 2 : length(x) 
        push!( X_hist, x[j] ) 
    end 
    t_hist = [ t_hist ; t_hist[end] .+ t[2:end] ]

    # set up next iter 
    xk = x[end]

# end 
X_hist = mapreduce( permutedims, vcat, X_hist ) 

## ============================================ ##
# propagate each segment 

X, t = prop_2Body_tof_Nseg( [r1; vz], dv_vec, N, tof_N, mu )  

## ============================================ ##
# compute miss distance 

dx_miss = miss_distance_prop2Body( x0, dv_vec, N, r2, tof_N, mu ) 

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
lines!( rv_lambert[:,1], rv_lambert[:,2], rv_lambert[:,3]; linewidth = 2, label = "lambert" ) 
scatter!( rv_lambert[1,1], rv_lambert[1,2], rv_lambert[1,3]; marker = :circle, markersize = 10, color = :black ) 
text!( rv_lambert[1,1], rv_lambert[1,2], rv_lambert[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 


# plot N segment solution 
lines!( X[:,1], X[:,2], X[:,3]; linewidth = 2, label = "lambert" ) 
scatter!( X[1,1], X[1,2], X[1,3]; marker = :circle, markersize = 10, color = :black ) 
# text!( X[1,1], X[1,2], X[1,3]; text = "P IC", color = :gray, offset = text_offset, align = (:center, :bottom) ) 

# target 
scatter!( r2[1], r2[2], r2[3]; marker = :circle, markersize = 10, color = :black )
text!( r2[1], r2[2], r2[3]; text = "target", color = :gray, offset = text_offset, align = ( :center, :bottom ) )  

display(fig) 
