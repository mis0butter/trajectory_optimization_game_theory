using LinearAlgebra



## ============================================ ##

# Function solves Kepler's equation M = E-e*sin(E)
# Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
# Output  eccentric anomaly E [rad]. 

export kepler_H 
function kepler_H( 
    M, 
    e, 
    eps = 1e-10, 
)

    H_k = M 
    H_kp1 = H_k + ( M - e*sinh(H_k) + H_k ) / ( e*cosh(H_k) - 1 ) 

    while abs(H_kp1 - H_k) > eps  
        H_k   = H_kp1 
        H_kp1 = H_k + ( M - e*sinh(H_k) + H_k ) / ( e*cosh(H_k) - 1 ) 
    end
    H = H_kp1 

    return H  
end 

## ============================================ ##

# Function solves Kepler's equation M = E-e*sin(E)
# Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
# Output  eccentric anomaly E [rad]. 

export kepler_E 
function kepler_E( 
    M, 
    e, 
    eps = 1e-10, 
)

    E_k  = M 
    E_kp1 = E_k - ( E_k - e*sin(E_k)- M ) / ( 1 - e*cos(E_k) ) 

    while ( abs(E_kp1-E_k) > eps )
        E_k   = E_kp1 
        E_kp1 = E_k - (E_k - e*sin(E_k) - M) / (1 - e*cos(E_k)) 
    end 
    E = E_kp1 

    return E 
end 

## ============================================ ##

export kepler_prop_tof 
function kepler_prop_tof( 
    rv_0, 
    tof, 
    mu = 1.0, 
)  

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
        
        # set target rv 
        oe_f    = copy(oe_0) 
        oe_f[6] = nu_f 
        rv_f    = kep2cart( oe_f, mu ) 

    else 

        println( "hyperbolic" ) 

        # initial hyperbolic anomaly 
        H_0 = acosh( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

        # mean motion 
        n = sqrt( complex(mu / (-a)^3) ) 
        
        # mean anomaly - propagate! 
        dM  = n * tof 
        M_0 = e * sinh(H_0) - H_0 
        M_f = M_0 + dM 

        # find final hyperbolic anomaly 
        H_f  = kepler_H( M_f, e ) 
        nu_f = 2*atan( sqrt( (e+1)/(e-1) ) * tanh(H_f/2) ) 

        # compute target rv 
        r_f = a*( 1 - e^2 ) / ( 1 + e * cos(nu_f) ) 
        v_f = sqrt( mu * ( 2/r_f - 1/a ) ) 
        rv_f = [r_f ; v_f] 

    end 

    return rv_f 
end 