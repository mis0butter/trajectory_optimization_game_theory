using trajectory_optimization_game_theory 
using ForwardDiff 

## ============================================ ##

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
rv_0 = [ r_0 ; v_0 ] 
rv_f = [ r_f ; v_f ] 
tof_N = tof / N 

miss_kepler = miss_distance_prop_kepler( 
    rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

obj( Δv_vec ) = miss_distance_prop_kepler( 
        rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

ForwardDiff.jacobian( obj, Δv_vec ) 
