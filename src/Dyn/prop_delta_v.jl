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
