

#============================================================

prop_2Body_dt_Nseg_range:

Description: Propagates an initial state through a vector of N trajectory segments

Inputs: 
    x0      - Initial state vector of form [r̄; v̄]
    Δv      - Matrix of size (N, 3) where each row is the velocity vector at a segment of the trajectory
    dt_N   - Change in time for each segment 
    N       - vector of segments of the trajectory
    mu      - Gravitational parameter (default = 1.0) 

Outputs:
    X_hist  - Matrix of size (N, 6) where each row is the state at each segment of the trajectory
    t       - Change in time between initial and final states

============================================================#

export prop_2Body_dt_Nseg 
function prop_2Body_dt_Nseg(
    x0, 
    Δv_vec, 
    N, 
    dt_N, 
    mu = 1.0 
    ) 
    
    # Creating Iteration Variables
    xk = copy(x0)
    Δt = N * dt_N 

    # Propagating Through Each Segment 
    X_hist = [ xk ] 
    t_hist = [ 0 ] 
    for i = 1 : N 

        # apply dv 
        xkdv = apply_dv( xk, Δv_vec[i,:] ) 

        # propagate and save 
        t, x = propagate_2Body( xkdv, dt_N, mu ) 
        for j = 2 : length(x) 
            push!( X_hist, x[j] ) 
        end 
        t_hist = [ t_hist ; t_hist[end] .+ t[2:end] ]

        # set up next iter 
        xk = x[end]

    end 
    X_hist = mapreduce( permutedims, vcat, X_hist ) 

    return X_hist, t_hist 
end


#============================================================

APPLY_DV:

Description: Helper function for prop_stateUV_Nseg() and prop_stateUV_Nseg_range() that adds Δv to a state vector's velocity

Inputs:
    1. x  - state vector
    2. Δv - velocity vector

Outputs:
    1. x - updated state vector

============================================================#

export apply_dv 
function apply_dv(x, Δv)

    # Applying Δv
    r = x[1:3]
    v = x[4:6] + Δv
    x = vcat(r, v)

    return x
end 


#============================================================
MISS_DISTANCE: 

Description: Calculates miss distance between trajectories

Inputs: 
    x0      - initial state vector of form [r; v]
    Δv_vec  - [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N       - number of segments 
    xf      - target state vector of form [r; v] 
    dt_N   - tof for each segment 
    mu      - Graviational parameter (default = 1.0) 

Outputs: 
    Δxf     - Miss distance
============================================================#

export miss_distance 
function miss_distance( x0, Δv_vec, N, xf, dt_N = 1.0, mu = 1.0 )

    # Propagating To Final State
    X_hist, Δt = prop_2Body_dt_Nseg( x0, Δv_vec, N, dt_N, mu ) 

    # extract final state 
    xf_prop = X_hist[end,:] 

    # Finding Miss Distance
    Δxf = abs.(xf_prop[1:3] - xf[1:3])

    return Δxf
end

