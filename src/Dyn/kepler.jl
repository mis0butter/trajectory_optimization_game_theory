using LinearAlgebra

## ============================================ ##

# Function solves Kepler's equation M = E-e*sin(E)
# Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
# Output  eccentric anomaly E [rad]. 

export keplter_E 
function kepler_E( 
    M, 
    e, 
    eps = 1e-10 
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
