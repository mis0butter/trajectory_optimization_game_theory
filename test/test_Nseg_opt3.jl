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
Δv_sol = min_Δv_dist( rv_0, rv_f, tof, N, mu ) 

# ----------------------- #

# create fig 
fig = plot_axes3d( )
fig = plot_orbit( rv_P, fig ) 
fig = plot_orbit( rv_E, fig ) 
fig = plot_prop_Δv( rv_0, Δv_sol, N, tof / N, mu, fig ) 





