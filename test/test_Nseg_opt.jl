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

x0_P = x_P[1] ;     r_0  = x0_P[1:3] ;  v_0 = x0_P[4:6] 
x0_E = x_E[1] ;     
xf_E = x_E[end] ;   r_f  = xf_E[1:3] ;  v_f = xf_E[4:6] 

x_P = vv2m(x_P) ;   
x_E = vv2m(x_E) ;   

## ============================================ ##
# test lambert soln 

# lambert solution 
dm      = "pro" 
tof     = t[end] / 2 
rv_prop_lb, Δv_vec  = prop_lambert_soln( x0_P, xf_E, tof, dm, mu )

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( rv_prop_lb, fig ) 
fig = plot_vector3d( [ x0_P[1:3] ], 500 * [ Δv_vec ], fig, :black, r/100 ) 


## ============================================ ##
# break up into 2 segments, see what happens 

rv_0_lb = rv_prop_lb[1,:] 
r_0_lb  = rv_0_lb[1:3] ;    v_0_lb  = rv_0_lb[4:6]  

tof_N = tof / 2 

# propagate first segment 
r_1         = r_0 
v_1         = v_f + Δv_vec / 2  
rv_1        = [ r_1 ; v_1 ] 
t, rv_prop_1 = propagate_2Body( rv_1, tof_N, mu ) 

# propagate second segment 
rv_2 = rv_prop_1[end] 
rv_2 = [ rv_2[1:3] ; rv_2[4:6] + Δv_vec / 2 ] 
t, rv_prop_2 = propagate_2Body( rv_2, tof_N, mu ) 

rv_prop_1 = vv2m( rv_prop_1 )
rv_prop_2 = vv2m( rv_prop_2 )  

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( rv_prop_1, fig ) 
fig = plot_orbit( rv_prop_2, fig ) 


