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

## ============================================ ##

"Calculate miss distance between trajectories using 2 body propagation" 
function miss_distance_prop2Body( 
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    rv_f,           # target state vector of form [r; v] 
    tof_N = 1.0,    # tof for each segment 
    mu = 1.0,       # gravitational parameter 
)

    # Propagating To Final State
    t, rv_hist = prop_2Body_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

    # extract final state 
    rv_f_prop = rv_hist[end,:] 

    # Finding Miss Distance
    Δrv_f = norm(rv_f_prop[1:3] - rv_f[1:3])

    return Δrv_f
end

export miss_distance_prop2Body 

## ============================================ ##

"Sum of Δv vector norms"
function sum_Δv_flat( tof_N_Δv_vec_flat, N )

    Δv_vec_flat = tof_N_Δv_vec_flat[2:end] 
    Δv_vec      = reshape( Δv_vec_flat, N, 3 ) 
    Δv          = [norm(Δv_vec[i, :]) for i in 1:N]
    fuel_norm   = sum(Δv)

    return fuel_norm
end

export sum_Δv_flat 

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

"Calculate miss distance between trajectories using kepler propagation" 
function miss_distance_prop_kepler( 
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    rv_f,           # target state vector of form [r; v] 
    tof_N = 1.0,    # tof for each segment 
    mu = 1.0,       # gravitational parameter 
)

    # Propagating To Final State
    t, rv_hist = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

    # extract final state 
    rv_f_prop = rv_hist[end,:] 

    # Finding Miss Distance
    Δrv_f = norm( rv_f_prop[1:3] - rv_f[1:3] ) 

    if isnan(Δrv_f) 
        println("Δrv_f is nan")
    end 

    return Δrv_f
end 

export miss_distance_prop_kepler 
# miss_kepler = miss_distance_prop_kepler( 
    # rv_0, Δv_vec, N, rv_f, tof_N, mu )

# ## ============================================ ##

"Compute miss distance between trajectories using kepler propagation and tof_N as part of the objective function with decision variable as [ tof_n ; Δv_flat ]"
function miss_tof_Δv_flat( 
    rv_0,               # initial state vector of form [r; v] 
    tof_N_Δv_flat,      # [N*3+1,1] vector of tof_N and Δv_flat 
    N,                  # number of segments 
    rv_f,               # target state vector of form [r; v] 
    mu = 1.0            # gravitational parameter 
) 

    # get tof and Δv 
    tof_N   = tof_N_Δv_flat[1] 
    Δv_flat = tof_N_Δv_flat[2:end] 
    
    Δv_vec  = reshape( Δv_flat, N, 3 ) 
    miss    = miss_distance_prop_kepler( rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
    return miss 
end 

export miss_tof_Δv_flat 

## ============================================ ##

"Compute sum of miss distance and magnitude of state vector"
function miss_mag_tof_Δv_flat( 
    rv_0,               # initial state vector of form [r; v] 
    tof_N_Δv_flat,      # [N*3+1,1] vector of tof_N and Δv_flat 
    N,                  # number of segments 
    rv_f,               # target state vector of form [r; v] 
    mu = 1.0            # gravitational parameter 
) 

    # get tof and Δv 
    tof_N       = tof_N_Δv_flat[1] 
    Δv_vec_flat = tof_N_Δv_flat[2:end] 
    
    Δv_vec = reshape( Δv_vec_flat, N, 3 ) 
    miss   = miss_distance_prop_kepler( rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

    # state magnitude 
    state_mag = norm( tof_N_Δv_flat ) 
    
    return miss + state_mag 
end 

export miss_mag_tof_Δv_flat 

