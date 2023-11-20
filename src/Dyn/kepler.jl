using LinearAlgebra

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

    En  = M 
    Ens = En - ( En - e*sin(En)- M ) / ( 1 - e*cos(En) ) 

    while ( abs(Ens-En) > eps )
        En  = Ens 
        Ens = En - (En - e*sin(En) - M) / (1 - e*cos(En)) 
    end 
    E = Ens 

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

    # eccentric anomaly 
    E_dv = acos( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

    # mean motion 
    n = sqrt( mu / a^3 ) 

    # mean anomaly 
    dM  = n * tof 
    M_0 = E_dv - e * sin(E_dv) 
    M_f = M_0 + dM 

    # eccentric anomaly 
    E_f = kepler_E( M_f, e ) 

    # get back true anomaly 
    nu_f = acos( ( cos(E_f) - e ) / ( 1 - e*cos(E_f) ) ) 

    # set target rv 
    oe_f    = copy(oe_0) 
    oe_f[6] = nu_f 
    rv_f    = kep2cart( oe_f, mu ) 

    return rv_f 
end 