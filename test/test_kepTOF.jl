using trajectory_optimization_game_theory
using LinearAlgebra 

## ============================================ ## 
# IC, target state, and lambert solve 

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 

# plot 
fig = plot_orbit( rv_lambert ) 

## ============================================ ##
# check Kepler TOF eqns --> given rv_0 and rv_f, TOF match 

tof_N = tof / N 
rv_0  = [ r_0; v_0 ]  
rv_k  = rv_0 

# apply delta v 
i = 1 
rv_0 = apply_Δv( rv_k, Δv_vec[i,:] ) 

# propagate using dynamics integration 
t, rv = propagate_2Body( rv_0, tof_N, mu ) 

# set target rv 
rv_f = rv[end] 

# get OE elements 
oe_0 = cart2kep( rv_0, mu ) 
oe_f = cart2kep( rv_f, mu ) 
a    = oe_0[1] 
e    = oe_0[2] 

# true anomaly 
nu_0 = oe_0[6] 
nu_f = oe_f[6] 
dnu  = nu_f - nu_0  

# eccentric anomaly 
E_dv = acos( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 
E_f  = acos( ( e + cos(nu_f) ) / ( 1 + e * cos(nu_f) ) ) 
dE   = E_f - E_dv  

# compute TOF 
TOF  = sqrt( a^3 / mu ) * ( E_f - e * sin(E_f) - E_dv + e * sin(E_dv) ) 

# check TOF and propagation t are same  
# println( "TOF = ", TOF ) 
# println( "t   = ", t[end] ) 
println( "TOF - t = ", TOF - t[end] ) 

## ============================================ ##
# check Kepler TOF eqns --> given rv_0 and TOF, rv_f match 

tof_N = tof / N 

# test for hyperbolic orbit 
rv_0  = [ r_0; v_0*1.3 ]  

# check that rv_0 is hyperbolic 
rv_k  = rv_0 
oe_k = cart2kep( rv_k, mu ) 
println( "e = ", oe_k[2] ) 

# propagate using dynamics integration 
t, rv = propagate_2Body( rv_0, tof_N, mu ) 

# convert to OEs 
oe_0 = cart2kep( rv_0, mu ) 
a    = oe_0[1] 
e    = oe_0[2] 
nu_0 = oe_0[6] 

if e < 1.0 

    # initial eccentric anomaly 
    E_0 = acos( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

    # mean motion 
    n = sqrt( mu / a^3 ) 

    # mean anomaly - propagate! 
    dM  = n * tof 
    M_0 = E_0 - e * sin(E_0) 
    M_f = M_0 + dM 
    
    # find final eccentric anomaly 
    E_f  = kepler_E( M_f, e ) 
    nu_f = acos( ( cos(E_f) - e ) / ( 1 - e*cos(E_f) ) ) 

else 

    println( "hyperbolic" ) 

    # initial hyperbolic anomaly 
    H_0 = acosh( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

    # mean motion 
    n = sqrt( mu / (abs(a))^3 ) 
    
    # mean anomaly - propagate! 
    dM  = n * tof 
    M_0 = e * sinh(H_0) - H_0 
    M_f = M_0 + dM 

    # find final hyperbolic anomaly 
    H_f  = kepler_H( M_f, e ) 
    nu_f = 2*atan( sqrt( (e+1)/(e-1) ) * tanh(H_f/2) ) 

end 

# set target rv 
oe_f    = copy(oe_0) 
oe_f[6] = nu_f 
rv_f    = kep2cart( oe_f, mu ) 

# check with function 
rv_f = prop_kepler_tof( rv_0, tof, mu )  

# println( "prop 2Body: rv[end] = ", rv[end] )
# println( "kep prop TOF: rv_f = ", rv_f )  
# println( "kep prop TOF check: rv_f_check = ", rv_f_check )  
println( "err norm = ", norm( rv[end] - rv_f ) ) 

## ============================================ ##
# use Kepler TOF equations to propagate each segment 

tof_N = tof / N 
rv_0  = [ r_0; v_0 ]  
rv_k  = rv_0 

rv_hist = [ apply_Δv( rv_k, Δv_vec[1,:] ) ] 
t_hist  = [ 0 ] 
for i = 1 : N 

    # apply delta v 
    rv_dv = apply_Δv( rv_k, Δv_vec[i,:] ) 

    # propagate using kepler TOF 
    rv_k = prop_kepler_tof( rv_dv, tof_N, mu ) 
    t, x = propagate_2Body( rv_dv, tof_N, mu ) 

    println( "prop err norm = ", norm( rv_k - x[end] ) ) 
    
    push!( rv_hist, rv_k ) 
    push!( t_hist, t_hist[end] + tof ) 

end 
rv_hist = mapreduce( permutedims, vcat, rv_hist ) 

## ============================================ ##
# test N segment propagation

t_2Body, rv_2Body = prop_2Body_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu )  
t_kep, rv_kep     = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

println( "final err norm = ", norm( rv_2Body[end,:] - rv_kep[end,:] ) )

plot_orbit( rv_2Body ) 

## ============================================ ##
# test miss distance 

rv_f = [r_f ; v_f]

miss_2Body = miss_distance_prop2Body( 
    rv_0, Δv_vec, N, rv_f, tof_N, mu)

miss_kepler = miss_distance_prop_kepler( 
rv_0, Δv_vec, N, rv_f, tof_N, mu)

## ============================================ ##
# now try ... minimizing delta v 

