using LinearAlgebra 

## ============================================ ##

"Create dummy IC for lambert transfer and then breaking up into smaller Δv. for testing purposes only!!!"
function lambert_IC() 

    # define IC, target state, and lambert solve 
    r_0 = [ 20.0e6, 20.0e6, 0]  # [m] 
    r_f = [-20.0e6, 10.0e6, 0]  # [m] 
    tof = 1.0 * 86400 
    mu  = 398600.4418e9         # [m^3/s^2] 
    dm  = "pro" 

    # solve lambert orbit 
    v_0, v_f = lambertbattin( r_0, r_f, mu, dm, tof ) 
    rv_0     = [ r_0; v_0 ] 

    # propagate lambert orbit 
    t_lambert, rv_lambert = propagate_2Body( rv_0, tof, mu, 1.0 ) 
    rv_lambert = mapreduce( permutedims, vcat, rv_lambert ) 

    # initial position has all z velocity 
    v_z = [ 0; 0; norm(v_0) ]

    # compute delta v vec for lambert solution 
    dv = v_0 - v_z   

    # break up delta v into smaller segments 
    N = 20 
    Δv_vec = [] 
    for i = 1 : N 
        push!( Δv_vec, dv/N ) 
    end 
    Δv_vec = mapreduce( permutedims, vcat, Δv_vec ) 

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
