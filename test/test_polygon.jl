using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 

## ============================================ ##
# init params 

mu = 398600.4415
r  = 6378.0
kep0_P = [ r+400.0, 0.1, -20*pi/180, 10.0*pi/180, 20.0*pi/180, 30.0*pi/180 ]
rv_0_P = kep2cart(kep0_P, mu) 
kep0_E = [ r+450.0, 0.2, 10.6*pi/180, 40.0*pi/180, 0.0, 120.0*pi/180 ]
rv_0_E = kep2cart(kep0_E, mu) 


# tof for pursuer to catch up to evader 
tof = 2000 

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

# plot 
fig = plot_axes3d( )
fig = plot_orbit( rv_E, fig ) 

## ============================================ ##
# ... it's just vector addition to find the vertices of the polygon 

# let's propagate the evader SC 
rv_f_E = rv_E[end,:] 

# center of polygon 
r_f = rv_f_E[1:3]  ; r̄_f = r_f / norm(r_f)  
v_f = rv_f_E[4:6]  ; v̄_f = v_f / norm(v_f) 

# define vector normal to orbit plane 
n̄_f = cross( r̄_f, v̄_f )  ; n̄_f = n̄_f / norm(n̄_f) 

# ok, let's plot this so that it all looks right 
fig = plot_vector3d( [ r_f ] , [ r̄_f * r ] , fig ) 
fig = plot_vector3d( [ r_f ] , [ n̄_f * r ] , fig ) 
fig = plot_vector3d( [ r_f ] , [ v̄_f * r ] , fig ) 


