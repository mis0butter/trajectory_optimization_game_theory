using LinearAlgebra 

## ============================================ ##

"Convert true anomaly to eccentric anomaly"
function nu2E( 
    nu,             # true anomaly [rad] 
    e,              # eccentricity 
) 


    # elliptic 
    if e < 1.0 
        E = 2 * atan( sqrt( (1-e)/(1+e) ) * tan(nu/2) ) 

    # hyperbolic 
    else 
        E = acosh( (e+cos(nu)) / (1+e*cos(nu)) ) 
        if nu < 0.0 || nu > pi 
            E = -E 
        end 
    end 

    if E < 0.0 
        E = 2*pi + E 
    end

    return E        # eccentric anomaly [rad] 
end 

export nu2E 

## ============================================ ##

"Convert eccentric anomaly to true anomaly" 
function E2nu( 
    E,              # eccentric anomaly [rad] 
    e,              # eccentricity 
) 

    # elliptic 
    if e < 1.0 
        nu = 2 * atan( sqrt( (1+e)/(1-e) ) * tan(E/2) ) 
    # hyperbolic 
    else 
        nu = acos( ( cosh(E) - e ) / ( 1 - e*cosh(E) ) ) 
        if E < 0.0 || E > pi 
            nu = -nu 
        end 
    end 

    if nu < 0.0 
        nu = 2*pi + nu 
    end 

    return nu       # true anomaly [rad] 
end 

export E2nu 

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
    E_k   = M 
    E_kp1 = E_k .- ( E_k .- e*sin(E_k) .- M ) ./ ( 1 .- e*cos(E_k) ) 

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
    # nu_f = nu_tof( e, n, nu_0, tof ) 
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

"Compute true anomaly from eccentric anomaly"
function nu_tof(  
    e,              # eccentricity 
    n,              # mean motion 
    nu_0,           # initial true anomaly 
    tof,            # time of flight 
) 

    # initial eccentric anomaly 
    E_0 = nu2E( nu_0, e ) 

    # mean anomaly - propagate! 
    dM  = n * tof 
    if e < 1.0 
        M_0 = E_0 - e * sin(E_0) 
    else 
        M_0 = e * sinh(E_0) - E_0 
    end 
    M_f = M_0 .+ dM 
    
    # find final eccentric anomaly 
    E_f  = kepler_E( M_f, e ) 
    nu_f = E2nu( E_f, e ) 

    return nu_f     # final true anomaly 
end 

export nu_tof 

# nu_f = elliptic_nu( e, n, nu_0, tof ) 

## ============================================ ##

"Compute true anomaly from eccentric anomaly for elliptic orbits"
function elliptic_nu(  
    e,              # eccentricity 
    n,              # mean motion 
    nu_0,           # initial true anomaly 
    tof,            # time of flight 
) 

    # initial eccentric anomaly 
    E_0 = nu2E( nu_0, e ) 
    # E_0 = acos( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

    # mean anomaly - propagate! 
    dM  = n * tof 
    M_0 = E_0 - e * sin(E_0) 
    M_f = M_0 .+ dM 
    
    # find final eccentric anomaly 
    E_f  = kepler_E( M_f, e ) 
    nu_f = E2nu( E_f, e ) 
    # nu_f = acos( ( cos(E_f) - e ) / ( 1 - e*cos(E_f) ) ) 

    return nu_f     # final true anomaly 
end 

export elliptic_nu 

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
    H_0 = nu2E( nu_0, e ) 
    # H_0 = acosh( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 
    
    # mean anomaly - propagate! 
    dM  = n * tof 
    M_0 = e * sinh(H_0) - H_0 
    M_f = M_0 + dM 

    # find final hyperbolic anomaly 
    H_f  = kepler_H( M_f, e ) 
    nu_f = E2nu( H_f, e )
    # nu_f = 2*atan( sqrt( (e+1)/(e-1) ) * tanh(H_f/2) ) 

    return nu_f     # final true anomaly 
end 

export hyperbolic_nu 

# nu_f = hyperbolic_nu( e, n, nu_0, tof ) 