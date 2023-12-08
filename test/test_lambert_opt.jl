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

x_P = vv2m(x_P) 
x_E = vv2m(x_E) 

## ============================================ ##
# plot 

r   = 6378.0
xyz = [ zeros(3) for i in 1:3 ] 
uvw = 0.5 * r .* [ [1,0,0] , [0,1,0] , [0,0,1] ] 

fig = plot_vector3d( [ xyz[1] ] , [ uvw[1] ], nothing, :red, r/100 ) 
fig = plot_vector3d( [ xyz[2] ] , [ uvw[2] ], fig, :blue, r/100 ) 
fig = plot_vector3d( [ xyz[3] ] , [ uvw[3] ], fig, :green, r/100 ) 

fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 





















## ============================================ ##
# init 

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
# rv_0  = [ r_0 ; v_0 ]
rv_0  = [ r_0 ; 0*v_0 ] 
rv_f  = [ r_f ; v_f ] 

## ============================================ ##

# can I take gradient of lambert? 

# fn(tof) = lambertbattin( r1, r2, mu, dm, tof )
# dfn     = tof --> ForwardDiff.derivative( fn, tof ) 

fn(x) = seebatt(x)
fn(rand(1)[1]) 
dfn = x -> ForwardDiff.derivative( fn, x )  
dfn( rand(1)[1] )

fn(x) = seebattk(x)
fn(rand(1)[1]) 
dfn = x -> ForwardDiff.derivative( fn, x )  
dfn( rand(1)[1] )

fn(x) = lambertbattin( r1, r2, mu, dm, x )
fn(rand(1)[1]) 
dfn = x -> ForwardDiff.derivative( fn, x )  
dfn( rand(1)[1] )

## ============================================ ##


