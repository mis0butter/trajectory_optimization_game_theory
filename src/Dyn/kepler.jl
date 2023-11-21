using LinearAlgebra

## ============================================ ##

"Solve Kepler's equation M = e*sinh(H) - H" 
function kepler_H( 
    M,              # mean anomaly [rad] 
    e,              # eccentricity 
    eps = 1e-10,    # tolerance 
)

    # starting guess 
    H_k   = M 
    H_kp1 = H_k + ( M - e*sinh(H_k) + H_k ) / ( e*cosh(H_k) - 1 ) 

    # Newton-Raphson iteration 
    while abs(H_kp1 - H_k) > eps  
        H_k   = H_kp1 
        H_kp1 = H_k + ( M - e*sinh(H_k) + H_k ) / ( e*cosh(H_k) - 1 ) 
    end
    H = H_kp1 

    return H        # hyperbolic anomaly [rad] 
end 

export kepler_H 
# H = kepler_H( M, e )

## ============================================ ## 

"Solve Kepler's equation M = E-e*sin(E)" 
function kepler_E( 
    M,              # mean anomaly [rad] 
    e,              # eccentricity 
    eps = 1e-10,    # tolerance 
)

    # starting guess 
    E_k  = M 
    E_kp1 = E_k - ( E_k - e*sin(E_k)- M ) / ( 1 - e*cos(E_k) ) 

    # Newton-Raphson iteration 
    while abs(E_kp1-E_k) > eps 
        E_k   = E_kp1 
        E_kp1 = E_k - (E_k - e*sin(E_k) - M) / (1 - e*cos(E_k)) 
    end 
    E = E_kp1 

    return E        # eccentric anomaly [rad] 
end 

export kepler_E 
# E = kepler_E( M, e )

## ============================================ ## 

# Resource: BMW 

"Propagate Keplerian orbit using TOF"
function prop_kepler_tof( 
    rv_0,           # initial position and velocity vectors 
    tof,            # time of flight 
    mu = 1.0,       # gravitational parameter 
)  

    # convert to OEs 
    oe_0 = cart2kep( rv_0, mu ) 
    a    = oe_0[1] 
    e    = oe_0[2] 
    nu_0 = oe_0[6] 
    
    # mean motion 
    n = sqrt( mu / (abs(a))^3 ) 

    # get true anomaly 
    if e < 1.0 
        nu_f = elliptic_nu( e, n, nu_0, tof ) 
    else 
        nu_f = hyperbolic_nu( e, n, nu_0, tof ) 
    end 

    # set target rv 
    oe_f    = copy(oe_0) 
    oe_f[6] = nu_f 
    rv_f    = kep2cart( oe_f, mu ) 

    return rv_f     # target rv 
end 

export prop_kepler_tof 
# rv_f = prop_keper_tof( rv_0, tof, mu ) 

## ============================================ ##

"Compute true anomaly from eccentric anomaly for elliptic orbits"
function elliptic_nu(  
    e,              # eccentricity 
    n,              # mean motion 
    nu_0,           # initial true anomaly 
    tof,            # time of flight 
) 

    # initial eccentric anomaly 
    E_0 = acos( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

    # mean anomaly - propagate! 
    dM  = n * tof 
    M_0 = E_0 - e * sin(E_0) 
    M_f = M_0 + dM 
    
    # find final eccentric anomaly 
    E_f  = kepler_E( M_f, e ) 
    nu_f = acos( ( cos(E_f) - e ) / ( 1 - e*cos(E_f) ) ) 

    return nu_f     # final true anomaly 
end 

# nu_f = elliptic_nu( e, n, nu_0, tof ) 

## ============================================ ##

"Compute true anomaly from hyperbolic anomaly for hyperbolic orbits"
function hyperbolic_nu(  
    e,              # eccentricity 
    n,              # mean motion 
    nu_0,           # initial true anomaly 
    tof,            # time of flight 
) 
    
    # initial hyperbolic anomaly 
    H_0 = acosh( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 
    
    # mean anomaly - propagate! 
    dM  = n * tof 
    M_0 = e * sinh(H_0) - H_0 
    M_f = M_0 + dM 

    # find final hyperbolic anomaly 
    H_f  = kepler_H( M_f, e ) 
    nu_f = 2*atan( sqrt( (e+1)/(e-1) ) * tanh(H_f/2) ) 

    return nu_f     # final true anomaly 
end

# nu_f = hyperbolic_nu( e, n, nu_0, tof ) 

## ============================================ ##

"Propagate an initial state through a vector of N trajectory segments using Kepler's equations" 
function prop_kepler_tof_Nseg(
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    tof_N,          # tof for each segment 
    mu = 1.0        # gravitational parameter 
) 
    
    # Creating Iteration Variables
    rv_k = copy(rv_0)

    # Propagating Through Each Segment 
    i = 1 
    rv_hist = [ apply_dv( rv_k, Δv_vec[i,:] ) ] 
    t_hist  = [ 0 ] 
    for i = 1 : N 

        # apply dv 
        rv_k_dv = apply_dv( rv_k, Δv_vec[i,:] ) 

        # propagate and save 
        rv_k = prop_kepler_tof( rv_k_dv, tof_N, mu ) 
        push!( rv_hist, rv_k ) 
        push!( t_hist, t_hist[end] + tof_N ) 

    end 
    rv_hist = mapreduce( permutedims, vcat, rv_hist ) 
 
    return  t_hist, 
            rv_hist 
end

export prop_kepler_tof_Nseg 
# t_kep, rv_kep = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

## ============================================ ##

"Calculate miss distance between trajectories using kepler propagation" 
function miss_distance_prop_kepler( 
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    rv_f,           # target state vector of form [r; v] 
    tof_N = 1.0,    # tof for each segment 
    mu = 1.0,       # gravitational parameter 
)

    # Propagating To Final State
    t, rv_hist = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

    # extract final state 
    rv_f_prop = rv_hist[end,:] 

    # Finding Miss Distance
    Δrv_f = abs.(rv_f_prop[1:3] - rv_f[1:3])

    return Δrv_f
end 

export miss_distance_prop_kepler 
# miss_kepler = miss_distance_prop_kepler( rv_0, Δv_vec, N, rv_f, tof_N, mu )

