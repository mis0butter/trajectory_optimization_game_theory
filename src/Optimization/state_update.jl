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
    x̄0::AbstractVector, 
    Δv̄::AbstractMatrix{T}, 
    UV::T, 
    N::Int
) where T<:Real # Iteration Variable

    # Creating Iteration Variables 
    x̄k = copy(x̄0)
    Δt = 0.0

    # Propagating to N 
    for i = 1:N 

        # Applying Δv
        x̄kdv = apply_dv(x̄k, Δv̄[i, :])

        # Propagating
        x̄k, δt = propKepUV(x̄kdv, UV) 

        # Updating
        Δt += δt

    end

    # Outputting
    return x̄k, Δt
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
    1. X̄ - Matrix of size (N, 6) where each row is the nondimensionalized state at each segment of the trajectory
    2. Δt - Change in time between initial and final states

============================================================#

function prop_stateUV_Nseg_range(
    x̄0::AbstractVector, 
    Δv̄::AbstractMatrix{T}, 
    UV::T, 
    N::UnitRange
) where T<:Real# Iteration Variable

    # Creating Iteration Variables
    x̄k = copy(x̄0)
    Δt = 0.0

    # Output Variable
    X̄ = zeros(last(N)+1, 6)
    X̄[1, :] = x̄k

    # Propagating to N
    for i in N
        # Applying Δv
        x̄kdv = apply_dv(x̄k, Δv̄[i, :])

        # Propagating
        x̄k, δt = propKepUV(x̄kdv, UV)

        # Updating
        X̄[i+1, :] = x̄k 
        Δt += δt
    end

    # Outputting
    return X̄, Δt
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
