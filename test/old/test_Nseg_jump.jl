using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 
using JuMP 
using Ipopt 

## ====================================================================
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
fig = plot_orbit( rv_P, fig ) 
fig = plot_orbit( rv_E, fig ) 

## ====================================================================
# break up into N segments, see what happens 

# define init and target vectors for pursuer 
rv_f = rv_E[end,:] 
rv_0 = rv_0_P 

N = 20 
dm = "pro"  # direction of motion

tof_N, Δv_vec = lambert_init_guess( rv_0, rv_f, tof, N, mu, dm ) 
x_0 = reshape( Δv_vec, N*3, 1 ) 

# JuMP-compatible miss distance function (uses sum of squares instead of norm)
# This avoids the norm() issue with JuMP variables
function miss_distance_prop2Body_jump(rv_0, x, N, rv_f, tof_N, mu)
    # Reshape x from (N*3,) vector to (N, 3) matrix
    Δv_vec = reshape(x, N, 3)
    
    # Propagate using 2-body (avoids kepler propagation which uses norm in cart2kep)
    t, rv_hist = prop_2Body_tof_Nseg(rv_0, Δv_vec, N, tof_N, mu)
    
    # Extract final state
    rv_f_prop = rv_hist[end, :]
    
    # Use sum of squares instead of norm (equivalent for minimization)
    Δrv_f_squared = sum((rv_f_prop[1:3] - rv_f[1:3]).^2)
    
    return Δrv_f_squared
end

# define objective function 
obj_fn(x) = sum_norm_Δv(x, N) + 
            miss_distance_prop2Body_jump(rv_0, x, N, rv_f, tof_N, mu)

# set model 
model = Model(Ipopt.Optimizer)
set_silent(model)

# define variables 
@variable(model, x[1:N*3])

# define objective 
@objective(model, Min, obj_fn(x)) 

# optimize 
optimize!(model)

# Extract solution
x_sol = value.(x)  # Get optimized variable values
Δv_sol = reshape(x_sol, N, 3)  # Reshape to (N, 3) matrix format

# Alternative methods (commented out):
# Δv_sol = min_Δv( rv_0, rv_f, tof, N, mu ) 
# Δv_sol = min_Δv_dist( rv_0, rv_f, tof, N, mu ) 
# Δv_sol = max_Δv_dist( rv_0, rv_f, tof, N, mu )

## ====================================================================

# create fig 
fig = plot_axes3d( )
fig = plot_orbit( rv_P, fig ) 
fig = plot_orbit( rv_E, fig ) 
fig = plot_prop_Δv( rv_0, Δv_sol, N, tof / N, mu, fig ) 


