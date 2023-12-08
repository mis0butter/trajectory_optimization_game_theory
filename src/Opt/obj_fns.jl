using LinearAlgebra 

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
    miss    = miss_distance_prop_kepler_Nseg( rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
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
    miss   = miss_distance_prop_kepler_Nseg( rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

    # state magnitude 
    state_mag = norm( tof_N_Δv_flat ) 
    
    return miss + state_mag 
end 

export miss_mag_tof_Δv_flat 

## ============================================ ##

"Calculate miss distance between trajectories using kepler propagation"
function miss_distance_prop_kepler( 
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    rv_f,           # target state vector of form [r; v] 
    tof = 1.0,      # tof for each segment 
    mu = 1.0,       # gravitational parameter 
)

    # add delta v to initial state 
    rv_Δv = rv_0 + [ zeros(3) ; Δv_vec]  

    # Propagating To Final State
    rv_f_prop = prop_kepler_tof( rv_Δv, tof, mu ) 

    # Finding Miss Distance
    Δrv_f = norm( rv_f_prop[1:3] - rv_f[1:3] ) 

    if isnan(Δrv_f) 
        println("Δrv_f is nan")
    end 

    return Δrv_f
end 

## ============================================ ##

"Calculate miss distance between trajectories with N segments using kepler propagation" 
function miss_distance_prop_kepler_Nseg( 
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

export miss_distance_prop_kepler_Nseg 
# miss_kepler = miss_distance_prop_kepler_Nseg( 
    # rv_0, Δv_vec, N, rv_f, tof_N, mu )