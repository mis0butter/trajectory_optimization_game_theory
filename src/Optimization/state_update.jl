#============================================================

prop_stateUV_Nseg:

Description: Propagates an initial state to the Nth segment of the trajectory

Inputs:
    1. x̄0 - Non-dimensionalized initial state vector of form [r̄; v̄]
    2. Δv̄ - Matrix of size (N, 3) where each row is the non-dimensionalized velocity vector at a segment of the trajectory
    3. UV - Kepler's Universal Variable
    4. N - Number of segments of the trajectory

Outputs:
    1. xk - Non-dimensionalized final state vector
    2. Δt - Change in time between states

============================================================#

function prop_stateUV_Nseg(
    x̄0, 
    Δv̄, 
    UV, 
    N
) 

    # Creating Iteration Variables 
    xk = copy(x̄0)
    Δt = 0.0

    # Propagating to N 
    for i = 1:N 

        # Applying Δv
        xkdv = apply_dv(xk, Δv̄[i, :])

        # Propagating
        xk, δt = propKepUV(xkdv, UV) 

        # Updating
        Δt += δt

    end

    # Outputting
    return xk, Δt
end

#============================================================

prop_stateUV_Nseg_RANGE:

Description: Propagates an initial state through a vector of N trajectory segments

Inputs: 
    1. x̄0 - Non-dimensionalized initial state vector of form [r̄; v̄]
    2. Δv̄ - Matrix of size (N, 3) where each row is the non-dimensionalized velocity vector at a segment of the trajectory
    3. UV - Kepler's Universal Variable
    4. N - vector of segments of the trajectory

Outputs:
    1. xk - Matrix of size (N, 6) where each row is the nondimensionalized state at each segment of the trajectory
    2. Δt - Change in time between initial and final states

============================================================#

function prop_stateUV_Nseg_range(
    x̄0, 
    Δv̄, 
    UV, 
    N
    ) 
    
    # Creating Iteration Variables
    xk = copy(x̄0)
    Δt = 0.0

    # Output Variable
    X = zeros(last(N)+1, 6)
    X[1, :] = xk

    # Propagating to N
    for i in N
        # Applying Δv
        xkdv = apply_dv(xk, Δv̄[i, :]) 

        # Propagating
        xk, δt = propKepUV(xkdv, UV)

        # Updating
        X[i+1, :] = xk
        Δt += δt
    end

    # Outputting
    return X, Δt
end

#============================================================

APPLY_DV:

Description: Helper function for prop_stateUV_Nseg() and prop_stateUV_Nseg_range() that adds Δv to a state vector's velocity

Inputs:
    1. x̄ - Non-Dimensionalized state vector
    2. Δv̄ - Non-Dimensionalized velocity vector

Outputs:
    1. x̄ - Updated state vector

============================================================#


function apply_dv(x̄, Δv̄)

    # Applying Δv
    r̄ = x̄[1:3]
    v̄ = x̄[4:6] + Δv̄
    x̄ = vcat(r̄, v̄)

    return x̄
end
