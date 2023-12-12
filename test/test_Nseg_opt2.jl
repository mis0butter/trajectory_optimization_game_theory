using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 

## ============================================ ##
# init 

mu = 398600.4415
r  = 6378.0
kep0_P = [r+400.0, 0.1, -20*pi/180, 10.0*pi/180, 20.0*pi/180, 30.0*pi/180]
kep0_E = [r+450.0, 0.2, 10.6*pi/180, 40.0*pi/180, 0.0, 120.0*pi/180]
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
rv_prop_lb, Δv  = prop_lambert_soln( x0_P, xf_E, tof, dm, mu )

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( rv_prop_lb, fig ) 
fig = plot_vector3d( [ x0_P[1:3] ], 500 * [ Δv ], fig ) 


## ============================================ ##
# break up into 2 segments, see what happens 

# set initial guess 
N       = 20 
tof_N   = tof / N / 4 

Δv_vec = [ Δv ]
for i = 1 : N-1 
    push!( Δv_vec, zeros(3) ) 
end 
Δv_vec      = vv2m( Δv_vec ) 
Δv_vec_flat = reshape( Δv_vec, N*3, 1 ) 
x_0         = 0.9 * [ tof_N ; Δv_vec_flat ]  

# define objective function 
obj_fn(x) = sum_norm_Δv( x, N ) 
obj_fn(x_0) 

# equality constraint 
c_fn(x) = miss_distance_prop_kepler_Nseg( rv_0, x[2:end], N, rv_f, x[1], mu ) 
c_fn(x_0) 

# inequality constraint ? 
Δv_max = 2.0 
h_fn(x) = constrain_Δv( x, N, Δv_max )
h_fn(x_0) 

# minimize constrained 
x_min  = min_aug_L( obj_fn, x_0, c_fn, h_fn ) 

# get solution 
Δv_sol    = reshape( x_min[2:end], N, 3 ) 
tof_N_sol = x_min[1] 

# ----------------------- #

# propagate 2 body 
t, rv_2Body = prop_2Body_tof_Nseg( rv_0, Δv_sol, N, tof_N_sol, mu ) 

# propagate kepler 
t, rv_kepler = prop_kepler_tof_Nseg( rv_0, Δv_sol, N, tof_N_sol, mu ) 

# plot 
fig = plot_axes3d( )
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( rv_2Body, fig ) 
# fig = plot_vector3d( [ x0_P[1:3] ], 500 * [ Δv ], fig ) 

# set up vector plotting 
nodes_N = rv_kepler[1:N, 1:3] 

xyz = copy(nodes_N) 
uvw = copy(2000 * Δv_sol)

fig = plot_vector3d( xyz, uvw, fig, 100 ) 




