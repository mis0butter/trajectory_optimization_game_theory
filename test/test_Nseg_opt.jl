using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 

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

x_P = vv2m(x_P) ;   r_0  = x0_P[1:3] 
x_E = vv2m(x_E) ;   r_f  = xf_E[1:3] 

## ============================================ ##
# test lambert soln 

# lambert solution 
dm      = "pro" 
tof     = t[end] / 2 
rv_prop = prop_lambert_soln( x0_P, xf_E, tof, dm, mu )

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( rv_prop, fig ) 
# fig = plot_scatter3d( xf_E[1], xf_E[2], xf_E[3], fig ) 
# fig = plot_scatter3d( rv_f_kepler[1], rv_f_kepler[2], rv_f_kepler[3], fig, :utriangle, :green ) 

