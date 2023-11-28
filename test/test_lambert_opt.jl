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
