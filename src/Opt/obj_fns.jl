using LinearAlgebra 

## ====================================================================

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

## ====================================================================

"Sum of Δv vector norms"
function sum_Δv_flat( tof_N_Δv_vec_flat, N )

    Δv_vec_flat = tof_N_Δv_vec_flat[2:end] 
    Δv_vec      = reshape( Δv_vec_flat, N, 3 ) 
    Δv          = [norm(Δv_vec[i, :]) for i in 1:N]
    fuel_norm   = sum(Δv)

    return fuel_norm
end

export sum_Δv_flat 

# ## ====================================================================

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

## ====================================================================

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

## ====================================================================

"Calculate miss distance between trajectories using kepler propagation"
function miss_distance_prop_kepler( 
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    rv_f,           # target state vector of form [r; v] 
    tof = 1.0,      # tof for each segment 
    mu  = 1.0,      # gravitational parameter 
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

## ====================================================================

"Calculate miss distance between trajectories with N segments using kepler propagation" 
function miss_distance_prop_kepler_Nseg( 
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    rv_f,           # target state vector of form [r; v] 
    tof_N = 1.0,    # tof for each segment 
    mu    = 1.0,    # gravitational parameter 
) 

    # check if Δv_vec is [N,3] 
    if size(Δv_vec) != (N,3) 
        Δv_vec = reshape( Δv_vec, N, 3 ) 
    end 

    # Propagating To Final State
    t, rv_hist = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

    # extract final state 
    rv_f_prop = rv_hist[end,:] 

    # Finding Miss Distance
    Δrv_f = norm( rv_f_prop[1:3] - rv_f[1:3] )

    # A non-finite value here is EXPECTED and handled: a large enough trial Δv makes the
    # Kepler solve blow up, and min_Δv_dist_solve's guarded objective substitutes a large
    # finite value so the line search backs off.  Printing it once per occurrence produced
    # millions of lines in a sweep, so the notice is gone; the guard is the real handling.

    return Δrv_f
end 

export miss_distance_prop_kepler_Nseg 
# miss_kepler = miss_distance_prop_kepler_Nseg( 
    # rv_0, Δv_vec, N, rv_f, tof_N, mu )

## ====================================================================
    
"""
    sum_norm_Δv(x, N)

Total ΔV: `Σᵢ ‖Δvᵢ‖`, in km/s.

Previously returned `Σᵢ ‖Δvᵢ‖²` — an energy-like quantity in (km/s)², not ΔV — despite the name
and docstring (defect D3b / D11). That mattered because this is summed against a terminal miss
distance in **km** inside `min_Δv_dist`: at Δv ≈ 0.01 km/s per segment the squared form is
~1e-3 while the miss term starts at 1–40 km, so fuel was outweighed roughly 1000:1 and was not
meaningfully optimized at all.

`ε` desingularizes the gradient: `‖v‖` is not differentiable at `v = 0`, and a zero Δv segment
is common here, so plain `norm` would hand ForwardDiff a NaN. With ε = 1e-12 km/s the
perturbation is ~1e-12 against segment norms of ~1e-2 — utterly negligible numerically, but it
keeps the objective smooth everywhere, which the gradient-based solvers in A6 require.
"""
function sum_norm_Δv( x, N; ε = 1e-12 )

    Δv_vec = reshape( x, N, 3 )

    sum_norm = zero(eltype(x))
    for i = 1 : N
        sum_norm += sqrt( sum( Δv_vec[i,:].^2 ) + ε^2 )
    end

    return sum_norm
end

export sum_norm_Δv 

## ====================================================================

"Inequality constraint: Δv vector norm <= Δv_max "
function constrain_Δv( 
    x,      # decision variable 
    N,      # number of segments 
    Δv_max  # maximum Δv magnitude 
) 

    Δv_vec = reshape( x, N, 3 ) 

    norm_Δv = [ ]
    for i = 1 : N 
        # sum_norm += norm( Δv_vec[i,:] ) 
        push!( norm_Δv, norm( Δv_vec[i,:] ) - Δv_max )
    end 

    return norm_Δv 
end 

export constrain_Δv 

    