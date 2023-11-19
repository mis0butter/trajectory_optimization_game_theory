using trajectory_optimization_game_theory
using LinearAlgebra 

## ============================================ ## 

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

fig = plot_orbit( rv_lambert ) 

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
# use Kepler TOF equations to propagate each segment 

i = 1 

rv = [ r1; v1 ]  
dt_N = tof / N 




