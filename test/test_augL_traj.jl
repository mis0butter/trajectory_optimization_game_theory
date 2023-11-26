using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences
using LinearAlgebra 

## ============================================ ##

# init 
r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
rv_0  = [ r_0 ; v_0 ] 
rv_f  = [ r_f ; v_f ] 

# reshape Δv_vec and add tof_N  
tof_N       = tof / N 
Δv_vec_flat = reshape( Δv_vec, N*3, 1 ) 

# set initial guess 
x_0 = [ tof_N ; Δv_vec_flat ]  

# objective function 
obj_fn(x) = norm( x ) 

# eq constraint function 
c_fn(x) = miss_tof_Δv_flat( rv_0, x, N, rv_f, mu ) 

# ineq constraint function 
h_fn(x) = -x[1] 

x_min = min_aug_L_eq_ineq( obj_fn, c_fn, h_fn, x_0, 100 )  




