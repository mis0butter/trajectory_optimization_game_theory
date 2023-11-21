using LinearAlgebra

## ============================================ ##

"Function solves Kepler's equation M = E-e*sin(E)" 
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

## ============================================ ##

"Function solves Kepler's equation M = E-e*sin(E)" 
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

## ============================================ ## 

# Resource: BMW 

"Propagate Keplerian orbit using TOF"
function kepler_prop_tof( 
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

    if e < 1.0 

        # initial eccentric anomaly 
        E_0 = acos( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

        # mean anomaly - propagate! 
        dM  = n * tof 
        M_0 = E_0 - e * sin(E_0) 
        M_f = M_0 + dM 
        
        # find final eccentric anomaly 
        E_f  = kepler_E( M_f, e ) 
        nu_f = acos( ( cos(E_f) - e ) / ( 1 - e*cos(E_f) ) ) 

    else 

        # initial hyperbolic anomaly 
        H_0 = acosh( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 
        
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

    return rv_f     # target rv 
end 

export kepler_prop_tof 
