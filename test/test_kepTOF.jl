using trajectory_optimization_game_theory
using LinearAlgebra 

## ============================================ ## 
# define IC, target state, and lambert solve 

r0      = [20.0e6, 20.0e6, 0]   # [m] 
rf      = [-20.0e6, 10.0e6, 0]  # [m] 
tof     = 1.0 * 86400 
mu      = 398600.4418e9         # [m^3/s^2] 
dm      = "pro" 
Dtsec   = tof 

## ============================================ ##
# solve and propagate lambert orbit 

v0, vf = lambertbattin( r0, rf, mu, dm, tof ) 

rv0 = [ r0; v0 ] 
t_lambert, rv_lambert = propagate_2Body( rv0, tof, mu, 1.0 ) 
rv_lambert = mapreduce( permutedims, vcat, rv_lambert ) 

## ============================================ ##
# plot 

fig = plot_orbit( rv_lambert ) 

## ============================================ ##
# break up delta v into smaller segments 

# initial position has all z velocity 
vz  = [ 0; 0; norm(v0) ]

# compute delta v vec for lambert solution 
dv = v0 - vz   

# break up delta v into smaller segments 
N = 10 
dv_vec = [] 
for i = 1 : N 
    push!( dv_vec, dv/N ) 
end 
dv_vec = mapreduce( permutedims, vcat, dv_vec ) 

## ============================================ ##
# check Kepler TOF eqns --> given rv0 and rvf, TOF match 

dt_N = tof / N 
rv0  = [ r0; v0 ]  
rvk  = rv0 

# apply delta v 
rv_dv = apply_dv( rvk, dv_vec[i,:] ) 

# propagate using dynamics integration 
t, rv = propagate_2Body( rv_dv, dt_N, mu ) 

# set target rv 
rvf = rv[end] 

# get OE elements 
oe_dv = cart2kep( rv_dv, mu ) 
oe_f  = cart2kep( rvf, mu ) 
a     = oe_dv[1] 
e     = oe_dv[2] 

# true anomaly 
nu_dv = oe_dv[6] 
nu_f  = oe_f[6] 
dnu   = nu_f - nu_dv  

# eccentric anomaly 
E_dv = acos( ( e + cos(nu_dv) ) / ( 1 + e * cos(nu_dv) ) ) 
E_f  = acos( ( e + cos(nu_f) ) / ( 1 + e * cos(nu_f) ) ) 
dE   = E_f - E_dv  

# compute TOF 
TOF  = sqrt( a^3 / mu ) * ( E_f - e * sin(E_f) - E_dv + e * sin(E_dv) ) 

# check TOF and propagation t are same  
println( "TOF = ", TOF ) 
println( "t   = ", t[end] ) 
println( "TOF - t = ", TOF - t[end] ) 

## ============================================ ##
# check Kepler TOF eqns --> given rv0 and TOF, rvf match 

dt_N = tof / N 
rv0  = [ r0; v0 ]  
rvk  = rv0 

# apply delta v 
rv_dv = apply_dv( rvk, dv_vec[i,:] ) 

# propagate using dynamics integration 
t, rv = propagate_2Body( rv_dv, dt_N, mu ) 

# convert to OEs 
oe_dv = cart2kep( rv_dv, mu ) 
a     = oe_dv[1] 
e     = oe_dv[2] 
nu_dv = oe_dv[6] 

# eccentric anomaly 
E_dv = acos( ( e + cos(nu_dv) ) / ( 1 + e * cos(nu_dv) ) ) 

# mean motion 
n = sqrt( mu / a^3 ) 

# mean anomaly 
dM   = n * dt_N 
M_dv = E_dv - e * sin(E_dv) 
M_f  = M_dv + dM 

# eccentric anomaly 
E_f = kepler_E( M_f, e ) 

# get back true anomaly 
nu_f = acos( ( cos(E_f) - e ) / ( 1 - e*cos(E_f) ) ) 

# set target rv 
oe_f = copy(oe_dv) 
oe_f[6] = nu_f 
rv_f = kep2cart( oe_f, mu ) 


## ============================================ ##
# use Kepler TOF equations to propagate each segment 

dt_N = tof / N 
rv0  = [ r0; v0 ]  
rvk  = rv0 

rv_hist = [ apply_dv( rvk, dv_vec[i,:] ) ] 
# for i = 1 : N 

    println( "i = ", i ) 

    # apply delta v 
    rv_dv = apply_dv( rvk, dv_vec[i,:] ) 

    # propagate using dynamics integration 
    t, x = propagate_2Body( rv_dv, dt_N, mu ) 

    # orbital elements 
    oek = cart2kep( rv_dv, mu ) 
    a   = oek[1]                # semimajor axis 
    e   = oek[2]                # eccentricity       
    n   = sqrt( mu / a^3 )      # mean motion 
    M   = n * dt_N              # mean anomaly 
    E   = kepler_E( M, e )      # eccentric anomaly  
    
    println( "e = ", e ) 
    println( "E = ", E ) 

    # compute nu from eccentric anomaly 
    # nu = 2 * atan( sqrt( (1+e)/(1-e) ) * tan(E/2) ) 
    nu = acos( ( cos(E) - e ) / ( 1 - e*cos(E) ) ) 
    oek[6] = nu 

    rvk = kep2cart( oek, mu ) 



    # push!( rv_hist, rvk ) 

# end 
rv_hist = mapreduce( permutedims, vcat, rv_hist ) 


