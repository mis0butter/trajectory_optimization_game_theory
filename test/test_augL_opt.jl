using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences
using LinearAlgebra 


## ============================================ ##
# init 

mu = 398600.4415
r  = 6378.0
kep0_P = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_E = [r+450.0, 0.0, 10.6*pi/180, 0.0, 0.0, 50.0*pi/180]
t = (0.0, 1*orbitPeriod(kep0_E, mu)) 
t_P, x_P = propagate_2Body(kep2cart(kep0_P, mu), t, mu, 1.0)
t_E, x_E = propagate_2Body(kep2cart(kep0_E, mu), t, mu, 1.0)

x0_P = x0_P_OG = x_P[1] 
x0_E = x0_E_OG = x_E[1] 
xf_E = xf_E_OG = x_E[end] 


## ============================================ ##
# lambert solution 

dm = "pro" 

tof = t[end] 

r1 = x0_P[1:3] 
r2 = xf_E[1:3] 
v1, v2 = lambertbattin(r1, r2, mu, dm, tof) 

x0_lambert = [r1 ; v1] 

Δv_vec = v2 - x0_E[4:6] 

# test function 
rv_f_kepler = prop_kepler_tof( x0_lambert, tof, mu ) 
t, rv_prop  = propagate_2Body( x0_lambert, tof, mu ) 
rv_prop = mapreduce( permutedims, vcat, rv_prop ) 

fig = plot_orbit( rv_prop ) 
fig = plot_scatter3d( xf_E[1], xf_E[2], xf_E[3], fig ) 
fig = plot_scatter3d( rv_f_kepler[1], rv_f_kepler[2], rv_f_kepler[3], fig, :utriangle, :green ) 

## ============================================ ##
# what is the true change in true anomaly? 

oe0_lambert = cart2kep( x0_lambert, mu ) 
oef_lambert = cart2kep( rv_prop[end,:], mu ) 
e  = oe0_lambert[2] 

# get true anomaly 
Δν = ( oef_lambert[6] - oe0_lambert[6] ) 
ν_0 = oe0_lambert[6] 
ν_f = oef_lambert[6]  

# get eccentric anomaly 
E_0 = acos( ( e + cosd(ν_0) ) / ( 1 + e * cosd(ν_0) ) ) 
E_f = acos( ( e + cosd(ν_f) ) / ( 1 + e * cosd(ν_f) ) ) 
ΔE  = E_f - E_0 

## ============================================ ##
# let's figure out what's going on 

rv_0 = x0_lambert 

# "Propagate Keplerian orbit using TOF"
# function prop_kepler_tof( 
#     rv_0,           # initial position and velocity vectors 
#     tof,            # time of flight 
#     mu = 1.0,       # gravitational parameter 
# )  

    # convert to OEs 
    oe_0 = cart2kep( rv_0, mu ) 
    a    = oe_0[1] 
    e    = oe_0[2] 
    nu_0 = oe_0[6] 
    
    # mean motion 
    n = sqrt( mu / (abs(a))^3 ) 

    # get true anomaly 
    #     nu_f = elliptic_nu( e, n, nu_0, tof ) 

        # initial eccentric anomaly 
        E_0 = acos( ( e + cos(nu_0) ) / ( 1 + e * cos(nu_0) ) ) 

        # mean anomaly - propagate! 
        dM  = n * tof 
        M_0 = E_0 - e * sin(E_0) 
        M_f = M_0 .+ dM 
        
        # find final eccentric anomaly 
        E_f  = kepler_E( M_f, e ) 
        nu_f = acos( ( cos(E_f) - e ) / ( 1 - e*cos(E_f) ) ) 

    # set target rv 
    oe_f    = copy(oe_0) 
    oe_f[6] = nu_f 
    rv_f    = kep2cart( oe_f, mu ) 

