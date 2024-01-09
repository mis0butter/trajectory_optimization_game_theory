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

x0_P = x_P[1] ;     r_0  = x0_P[1:3] ;  v_0 = x0_P[4:6] ;   
rv_0 = [ r_0 ; v_0 ] 
x0_E = x_E[1] ;     
xf_E = x_E[end] ;   r_f  = xf_E[1:3] ;  v_f = xf_E[4:6] 
rv_f = [ r_f ; v_f ] 

x_P = vv2m(x_P) ;   
x_E = vv2m(x_E) ;   

## ============================================ ##
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

tof_N = tof / 2

# propagate first segment 
r_1          = r_0 
v_1          = v_0 + Δv 
rv_1         = [ r_1 ; v_1 ] 
t, rv_prop_1 = propagate_2Body( rv_1, tof_N, mu ) 

# propagate second segment 
rv_2 = rv_prop_1[end] 
rv_2 = [ rv_2[1:3] ; rv_2[4:6] + 0 * Δv ] 
t, rv_prop_2 = propagate_2Body( rv_2, tof_N, mu ) 

rv_prop_1 = vv2m( rv_prop_1 )
rv_prop_2 = vv2m( rv_prop_2 )  

## ============================================ ##
# now vary segments 

N       = 2 
tof_N   = tof / N 
Δv_vec  = vv2m( [ Δv, zeros(3) ] ) 
t_kep, rv_kep = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

miss_kepler = miss_distance_prop_kepler_Nseg( 
    rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

    
## ============================================ ##
# can I use Optim? ... looks like YES 

x_0 = reshape( 1.1 * Δv_vec, N*3, 1 ) 

# optimization 
fn(x)       = miss_distance_prop_kepler_Nseg( rv_0, x, N, rv_f, tof_N, mu ) 
x_min       = min_optim( fn, x_0 ) 
fn(x_min) 

# test 
Δv_vec = reshape( x_min, N, 3 ) 
t, rv_2Body = prop_2Body_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

## ============================================ ##
# ok, now let's try minimizing the miss distance for N segments 

# init stuff 
N       = 2 
tof_N   = tof / N 

# set initial guess 
Δv_vec_flat = reshape( Δv_vec, N*3, 1 ) 
x_0         = 1.1 * [ tof_N ; Δv_vec_flat ]  

obj_fn(x) = sum( abs.( x[2:end] ) ) 

x_min = min_aug_L( obj_fn, x_0  ) 

# equality constraint 
c_fn(x) = miss_distance_prop_kepler_Nseg( rv_0, x[2:end], N, rv_f, x[1], mu ) 
c_fn(x_0) 

# equality-constrained 
x_min = min_aug_L( obj_fn, x_0, c_fn )

# test 
Δv_vec = reshape( x_min[2:end], N, 3 ) 
tof_N  = x_min[1] 
t, rv_2Body = prop_2Body_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( rv_2Body, fig ) 
# fig = plot_vector3d( [ x0_P[1:3] ], 500 * [ Δv ], fig ) 



