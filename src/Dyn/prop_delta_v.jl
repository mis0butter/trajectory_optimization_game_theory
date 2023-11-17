

#============================================================

prop_state_dt_Nseg_range:

Description: Propagates an initial state through a vector of N trajectory segments

Inputs: 
    x0 - Initial state vector of form [r̄; v̄]
    Δv - Matrix of size (N, 3) where each row is the non-dimensionalized velocity vector at a segment of the trajectory
    dt - Change in time for each segment 
    N  - vector of segments of the trajectory
    mu - Gravitational parameter (default = 1.0) 

Outputs:
    X  - Matrix of size (N, 6) where each row is the state at each segment of the trajectory
    Δt - Change in time between initial and final states

============================================================#

export prop_state_dt_Nseg_range 
function prop_state_dt_Nseg_range(
    x0, 
    Δv, 
    dt, 
    N, 
    mu = 1.0 
    ) 
    
    # Creating Iteration Variables
    xk = copy(x0)
    Δt = N * dt 

    # Output Variable
    X = zeros(last(N)+1, 6)
    X[1, :] = xk

    # propagate for N segments 
    for i in 1:N

        # apply Δv
        xkdv = apply_dv(xk, Δv[i, :]) 

        # propagate segment 
        # xk, δt = propKepUV(xkdv, UV)
        t, x = propagate_2Body( xkdv, dt, mu) 

        # update final state 
        xk = x[end] 
        X[i+1, :] = xk 

    end

    return X, Δt
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

