function prop_opt_Nseg( 
    rv_0, 
    rv_f, 
    tof, 
    dm, 
    N  = 10, 
    mu = 1.0 
) 

    # first compute lambert 
    _, Δv  = prop_lambert_soln( rv_0, rv_f, tof, dm, mu )

    # set initial guess 
    tof_N   = tof / N / 4 

    Δv_vec = [ Δv ]
    for i = 1 : N-1 
        push!( Δv_vec, zeros(3) ) 
    end 
    Δv_vec      = vv2m( Δv_vec ) 
    Δv_vec_flat = reshape( Δv_vec, N*3, 1 ) 
    x_0         = 0.9 * [ tof_N ; Δv_vec_flat ]  

    # define objective function 
    obj_fn(x) = sum_norm_Δv( x, N ) 
    obj_fn(x_0) 

    # equality constraint 
    c_fn(x) = miss_distance_prop_kepler_Nseg( rv_0, x[2:end], N, rv_f, x[1], mu ) 
    c_fn(x_0) 

    # inequality constraint ? 
    Δv_max = 2.0 
    h_fn(x) = constrain_Δv( x, N, Δv_max )
    h_fn(x_0) 

    # minimize constrained 
    x_min  = min_aug_L( obj_fn, x_0, c_fn, h_fn ) 

    # get solution 
    Δv_sol    = reshape( x_min[2:end], N, 3 ) 
    tof_N_sol = x_min[1] 

    return tof_N_sol, Δv_sol 
end 

## ============================================ ##

"Propagates an initial state through a vector of N trajectory segments using dynamics integration" 
function prop_2Body_tof_Nseg(
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    tof_N,          # tof for each segment 
    mu = 1.0        # gravitational parameter 
) 
    
    # Creating Iteration Variables
    rv_k = copy(rv_0)
    Δt   = N * tof_N 

    # Propagating Through Each Segment 
    rv_hist = [ rv_k ] 
    t_hist  = [ 0 ] 
    for i = 1 : N 

        # apply dv 
        rv_k_dv = apply_Δv( rv_k, Δv_vec[i,:] ) 

        # propagate and save 
        t, rv = propagate_2Body( rv_k_dv, tof_N, mu ) 
        for j = 2 : length(rv) 
            push!( rv_hist, rv[j] ) 
        end 
        t_hist = [ t_hist ; t_hist[end] .+ t[2:end] ]

        # set up next iter 
        rv_k = rv[end]

    end 
    rv_hist = mapreduce( permutedims, vcat, rv_hist ) 
 
    return t_hist, rv_hist 
end

export prop_2Body_tof_Nseg 

## ============================================ ## 

"Propagate an initial state through a vector of N trajectory segments using Kepler's equations" 
function prop_kepler_tof_Nseg(
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    tof_N,          # tof for each segment 
    mu = 1.0        # gravitational parameter 
) 
    
    # set up time and state hists 
    rv_hist = [ apply_Δv( rv_0, Δv_vec[1,:] ) ] 
    t_hist  = [ ] ; push!( t_hist, 0.0 ) 

    # propagate Through Each Segment 
    rv_k = copy(rv_0)
    for i = 1 : N 

        # apply dv and propagate 
        rv_k_dv = apply_Δv( rv_k, Δv_vec[i,:] ) 
        rv_k    = prop_kepler_tof( rv_k_dv, tof_N, mu ) 

        # save 
        push!( rv_hist, rv_k ) 
        push!( t_hist, t_hist[end] .+ tof_N ) 

    end 
    rv_hist = mapreduce( permutedims, vcat, rv_hist ) 
 
    return  t_hist, 
            rv_hist 
end

export prop_kepler_tof_Nseg 
# t_kep, rv_kep = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

## ============================================ ##

"Adds Δv to a state vector's velocity" 
function apply_Δv(
    rv,         # state vector 
    Δv,         # velocity vector 
)

    # Applying Δv
    r = rv[1:3]
    v = rv[4:6] + Δv
    rv = vcat(r, v)

    return rv
end 

export apply_Δv 

