using trajectory_optimization_game_theory
using LinearAlgebra 

## ============================================ ## 

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
dv  = v0 - vz   

# break up delta v into smaller segments 
N = 10 
dv_vec = [] 
for i = 1 : N 
    push!( dv_vec, dv/N ) 
end 
dv_vec = mapreduce( permutedims, vcat, dv_vec ) 

## ============================================ ##
# use Kepler TOF equations to propagate each segment 

dt_N = tof / N 
rv0  = [ r0; v0 ]  
rvk  = rv0 

rv_hist = [ rvk + [ zeros(3) ; dv_vec[i,:] ] ] 
for i = 1 : N 

    println( "i = ", i ) 

    # apply delta v 
    rvk[4:6] = rvk[4:6] + dv_vec[i,:] 

    # orbital elements 
    oek = cart2kep( rvk, mu ) 
    a   = oek[1]                # semimajor axis 
    e   = oek[2]                # eccentricity       
    n   = sqrt( mu / a^3 )      # mean motion 
    M   = n * dt_N              # mean anomaly 
    E   = kepler_E( M, e )    # eccentric anomaly  
    
    println( "e = ", e ) 
    println( "E = ", E ) 

    # compute nu from eccentric anomaly 
    # nu = 2 * atan( sqrt( (1+e)/(1-e) ) * tan(E/2) ) 
    nu = acos( ( cos(E) - e ) / ( 1 - e*cos(E) ) ) 
    oek[6] = nu 

    rvk = kep2cart( oek, mu ) 
    push!( rv_hist, rvk ) 

end 
rv_hist = mapreduce( permutedims, vcat, rv_hist ) 



## ============================================ ##

# Function solves Kepler's equation M = E-e*sin(E)
# Input - Mean anomaly M [rad] , Eccentricity e and Epsilon 
# Output  eccentric anomaly E [rad]. 
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



