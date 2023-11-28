using LinearAlgebra 

## ============================================ ##

"Create dummy IC for lambert transfer and then breaking up into smaller Δv. for testing purposes only!!!"
function lambert_IC() 

    # define IC, target state, and lambert solve 
    R   = 6378.0                # Earth radius [km] 
    r_0 = [ 20.0e6, 20.0e6, 0]  # [km] 
    r_f = [-20.0e6, 10.0e6, 0]  # [km] 
    # r_0 = [ R + 500, R + 500, 0 ]
    # r_f = [ R - 500, R - 500, 0 ] 
    tof = 1.0 * 86400 
    mu  = 398600.4418e9         # [km^3/s^2] 
    dm  = "pro" 

    # solve lambert orbit 
    v_0, v_f = lambertbattin( r_0, r_f, mu, dm, tof ) 

    # # # let's non-dimensionalize everything 
    # r_0, v_0 = nondim_rv( r_0, v_0, mu, R )
    # r_f, v_f = nondim_rv( r_f, v_f, mu, R )  
    # mu = 1.0 

    rv_0     = [ r_0; v_0 ] 

    # propagate lambert orbit 
    t_lambert, rv_lambert = propagate_2Body( rv_0, tof, mu, 1.0 ) 
    rv_lambert = mapreduce( permutedims, vcat, rv_lambert ) 

    # N segments 
    N = 20 
    Δv_vec = zeros(N, 3) 
    Δv_vec[1,:] = v_0 

    # # initial position has all z velocity 
    # v_z = [ 0; 0; norm(v_0) ]

    # # compute delta v vec for lambert solution 
    # dv = v_0 - v_z   
    # dv = v_0 

    # # break up delta v into smaller segments 
    # Δv_vec = [] 
    # for i = 1 : N 
    #     push!( Δv_vec, dv/N ) 
    # end 
    # Δv_vec = mapreduce( permutedims, vcat, Δv_vec ) 

    return  r_0,            # initial position 
            r_f,            # final position 
            v_0,            # initial velocity 
            v_f,            # final velocity 
            rv_lambert,     # lambert orbit 
            Δv_vec,         # delta v vector 
            tof,            # time of flight 
            N,              # number of segments 
            mu              # gravitational parameter  
end 

export lambert_IC 
# r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 

