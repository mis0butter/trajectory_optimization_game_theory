using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 

## ============================================ ##
# init 

mu = 398600.4415
r  = 6378.0
kep0_P = [ r+400.0, 0.1, -20*pi/180, 10.0*pi/180, 20.0*pi/180, 30.0*pi/180 ]
kep0_E = [ r+450.0, 0.2, 10.6*pi/180, 40.0*pi/180, 0.0, 40.0*pi/180 ]
t = (0.0, 1*orbitPeriod(kep0_P, mu)) 
t_P, x_P = propagate_2Body(kep2cart(kep0_P, mu), t, mu, 1.0)
t = (0.0, 1*orbitPeriod(kep0_E, mu)) 
t_E, x_E = propagate_2Body(kep2cart(kep0_E, mu), t, mu, 1.0)

x0_P = x_P[1] ;     r_0  = x0_P[1:3] ;  v_0 = x0_P[4:6] ;   
rv_0 = [ r_0 ; v_0 ] 
x0_E = x_E[1] ;     
xf_E = x_E[end] ;   r_f  = xf_E[1:3] ;  v_f = xf_E[4:6] 
rv_f = [ r_f ; v_f ] 

x_P = vv2m(x_P) ;   
x_E = vv2m(x_E) ;   

# ----------------------- # 
# test lambert soln 

# lambert solution 
dm      = "pro" 
tof     = t[end] / 2 
rv_prop_lb, Δv  = prop_lambert_soln( rv_0, rv_f, tof, dm, mu )

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( rv_prop_lb, fig ) 
fig = plot_vector3d( [ rv_0[1:3] ], 500 * [ Δv ], fig ) 


## ============================================ ##

# init params 
mu = 398600.4415
r  = 6378.0
kep0_P = [ r+400.0, 0.1, -20*pi/180, 10.0*pi/180, 20.0*pi/180, 30.0*pi/180 ]
rv_0_P = kep2cart(kep0_P, mu) 
kep0_E = [ r+450.0, 0.2, 10.6*pi/180, 40.0*pi/180, 0.0, 40.0*pi/180 ]
rv_0_E = kep2cart(kep0_E, mu) 

# tof for pursuer to catch up to evader 
tof = 2000 

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

# plot 
fig = plot_axes3d( )
fig = plot_orbit( rv_P, fig ) 
fig = plot_orbit( rv_E, fig ) 


## ============================================ ##
# break up into N segments, see what happens 

# define init and target vectors for pursuer 
rv_f = rv_E[end,:] 
rv_0 = rv_0_P 

N = 20 
Δv_sol = min_Δv( rv_0, rv_f, tof, N, mu ) 

# ----------------------- #

# create fig 
fig = plot_axes3d( )
fig = plot_orbit( rv_P, fig ) 
fig = plot_orbit( rv_E, fig ) 
fig = plot_prop_Δv( rv_0, Δv_sol, N, tof / N, mu, fig ) 





