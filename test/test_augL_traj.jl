using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences
using LinearAlgebra 

## ============================================ ##
# init 

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
# rv_0  = [ r_0 ; v_0 ]
rv_0  = [ r_0 ; 0*v_0 ] 
rv_f  = [ r_f ; v_f ] 

# reshape Δv_vec and add tof_N  
tof_N       = tof / N 
# tof_N       = tof 
Δv_vec_flat = reshape( Δv_vec, N*3, 1 ) 

# set initial guess 
x_0 = [ tof_N ; Δv_vec_flat ]  

# plot solution 
t_kep, rv_kep_fd = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N_fd, mu ) 
fig = plot_orbit( rv_kep_fd ) 
fig = plot_scatter3d( rv_kep_fd[:,1], rv_kep_fd[:,2], rv_kep_fd[:,3], fig ) 
fig = plot_scatter3d( r_f[1], r_f[2], r_f[3], fig ) 


## ============================================ ##
# test BFGS directly 

# define objective fn 
fn( tof_Δv ) = miss_mag_tof_Δv_flat( rv_0, tof_Δv, N, rv_f, mu ) 
fn( x_0 ) 

# create gradient fn  
dfn_fd  = tof_Δv -> ForwardDiff.gradient( fn, tof_Δv ) 

# compute gradient 
g_fd  = dfn_fd( x_0 ) 

# minimize 
x_min_fd   = min_bfgs( fn, dfn_fd, x_0 )  

# get soln 
tof_N_fd  = x_min_fd[1] 
Δv_sol_fd = reshape( x_min_fd[2:end], N, 3 ) 

# plot solution 
t_kep, rv_kep_fd = prop_kepler_tof_Nseg( rv_0, Δv_sol_fd, N, tof_N_fd, mu ) 
fig = plot_orbit( rv_kep_fd ) 
fig = plot_scatter3d( rv_kep_fd[:,1], rv_kep_fd[:,2], rv_kep_fd[:,3], fig ) 
fig = plot_scatter3d( r_f[1], r_f[2], r_f[3], fig ) 

## ============================================ ##
## ============================================ ##
# inequality-constrained augmented Lagrangian 

# define objective fn 
obj_fn( tof_Δv ) = miss_mag_tof_Δv_flat( rv_0, tof_Δv, N, rv_f, mu ) 
obj_fn( x_0 ) 

# ineq constraint function 
h_fn(x) = -x[1] 
N_h     = length( h_fn(x_0) ) 

# lagrange multipliers and penalty parameters 
λ_0 = zeros(N_h) 
p_0 = 10.0 * ones(N_h) 
γ   = 2.0 

# define tol 
tol = 1e-6 

# ----------------------- #

# step 0: initialize 
λ_k = copy( λ_0 )
p_k = copy( p_0 ) 
x_k = copy( x_0 ) 
h_k = h_fn( x_k )

k = 0 ; loop = true 
while loop 

    k += 1 

    # step 1: assign augmented Lagrangian fn 
    fn(x_k) = aug_L_fn( obj_fn, h_fn, x_k, λ_k, p_k ) 
    dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

    # step 2: minimize unconstrained problem  
    x_min = min_bfgs( fn, dfn, x_k )  

    # step 3 check convergence ... 
    dx = norm(x_min - x_k) 
    if dx < tol 
        loop = false 
    else 
        x_k = x_min 
    end 

    # update constraint values 
    h_k = h_fn( x_k )
    N_h = length( h_k ) 

    # ... and update parameters 
    if N_h == 1 

        # update λ 
        λ_k = max( λ_k + p_k * h_k , 0.0 ) 

        # update p 
        if h_k > 0 
            p_k *= γ 
        end 
        
    else 

        for i in eachindex(λ_k)

            # update λ
            λ_k[i] = max( λ_k[i] + p_k[i] * h_k[i] , 0.0 ) 
    
            # update p 
            if h_k[i] > 0 
                p_k[i] *= γ 
            end 
    
        end    

    end 

end 

println( "x min = ", x_k ) 


## ============================================ ##
## ============================================ ##
# augmented Lagrangian ... equality-constrained 

# objective function 
obj_fn(x) = sum_Δv_flat( x, N ) 
# obj_fn(x) = norm( x[2:end] )^2 
# obj_fn(x) = miss_mag_tof_Δv_flat( rv_0, x, N, rv_f, mu ) 

# eq constraint function 
c_fn(x) = miss_tof_Δv_flat( rv_0, x, N, rv_f, mu ) 
N_c     = length( c_fn(x_0) ) 

# ----------------------- #
# augmented Lagrangian method (equality-constrained) 

# step 0: initialize 
λ_k = copy( λ_0 )
p_k = copy( p_0 ) 
x_k = copy( x_0 ) 

k = 0 ; loop = true 
while loop 

    k += 1 

    # step 1: assign augmented Lagrangian fn 
    fn(x_k) = aug_L_fn( obj_fn, c_fn, x_k, λ_k, p_k ) 
    dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

    # step 2: minimize unconstrained problem  
    x_min = min_bfgs( fn, dfn, x_k )  
    
    # step 3 check convergence ... 
    dx = norm(x_min - x_k) 
    if dx < tol 
        loop = false 
    else 
        x_k = x_min 
    end 

    # step 3: check constraint function and update parameters 
    if norm( c_fn(x_k) ) > tol 
        λ_k += p_k .* c_fn(x_k) 
        p_k *= γ 
    else 
        loop = false 
    end 

end 









## ============================================ ##
## ============================================ ##

# objective function 
obj_fn(x) = sum_Δv_flat( x, N ) 
# obj_fn(x) = miss_mag_tof_Δv_flat( rv_0, x, N, rv_f, mu ) 

# eq constraint function 
c_fn(x) = miss_tof_Δv_flat( rv_0, x, N, rv_f, mu ) 

# ineq constraint function 
h_fn(x) = -x[1] 

# minimize 
# x_min = min_aug_L_eq( obj_fn, c_fn, x_0, 100 )  
x_min = min_aug_L_eq_ineq( obj_fn, c_fn, h_fn, x_0, 100 )  
# x_min = min_aug_L_ineq( obj_fn, h_fn, x_0, 10 ) 




