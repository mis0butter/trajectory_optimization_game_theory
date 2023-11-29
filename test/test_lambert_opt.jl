using trajectory_optimization_game_theory
using LinearAlgebra 
using ForwardDiff 
using GLMakie 

## ============================================ ##
# example test 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# obj fn - true min at [-1, 0, 0]
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# eq constraints 
c_fn(x) = x[3] - 1          # x[3] - 1 = 0  

# ineq constraints: h_fn formulated as <= 0 
h_fn(x) = [ -x[1] + 1 ;     # x[1] >= 1   * active 
            -x[2] - 1 ]     # x[2] >= -1  * inactive 

x_min = min_aug_L_eq_ineq( obj_fn, c_fn, h_fn, x_0 ) 

## ============================================ ##
# 
# intercept problem: 
# 
# minimize 
# 		delta_v_1 
# 
# subject to: 	
# 		tof = min( tof_pro, tof_retro ) 
# 		v0  = min( v0_pro, v0_retro ) 
# 
# 		tof_pro,   v0_pro   = lambert( r_init, r_target, 'pro' )
# 		tof_retro, v0_retro = lambert( r_init, r_target, 'pro' )
# 
# 		delta_v <= constraint 
# 
# state: 		
# 		tof, v0 

tof = 1.0 
v0  = [ 1.0, 0.0, 0.0 ] 



