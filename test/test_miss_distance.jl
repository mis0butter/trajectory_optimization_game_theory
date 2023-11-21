using trajectory_optimization_game_theory 

## ============================================ ##

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof = lambert_IC() 

## ============================================ ##
# test miss distance 

rv_f = [r_f ; v_f]

miss_2Body = miss_distance_prop2Body( 
    rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

miss_kepler = miss_distance_prop_kepler( 
    rv_0, Δv_vec, N, rv_f, tof_N, mu )

