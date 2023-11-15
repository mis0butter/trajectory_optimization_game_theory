#============================================================
EQUALITY_CONSTRAINTS:

Description: Creates an equality constraints vector based off of propagation

Inputs:
    1. Δτ - Kepler's Universal Variable
    2. Δv_vec - Matrix of size (N, 3) where each row is the velocity vector at a segment of the trajectory
    3. x̄₀ - Non-dimensional initial state vector
    4. x̄f₀ - Non-dimensional final state vector
    5. N - Number of segments of the trajectory

Outputs:
    1. h_vec - Equality constraints vector
============================================================#

function equality_constraints(
    Δτ, 
    Δv_vec, 
    x̄₀, 
    x̄f₀, 
    N
)

    h_vec = zeros(3)
    
    # Terminal Constraint
    Δr = calculate_miss(Δτ, Δv_vec, x̄₀, x̄f₀, N)
    h_vec[1:3] = Δr

    return h_vec
end

#============================================================
INEQUALITY_CONSTRAINTS:

Description: Creates inequality constraints vector based off of velocity

Inputs:
    1. Δτ     - Kepler's Universal Variable
    2. Δv_vec - Matrix of size (N, 3) where each row is the velocity vector at a segment of the trajectory
    3. x̄₀     - Non-dimensional initial state vector
    4. x̄f₀    - Non-dimensional final state vector
    5. N      - Number of segments of the trajectory

Outputs:
    1. ψ_vec - Inequality constraints vector
============================================================#

function inequality_constraints(
    Δτ, 
    Δv_vec, 
    x̄₀, 
    x̄f₀, 
    N
)

    ψ_vec = zeros(N+1)

    # Δτ > 0 Inequality
    ψ_vec[1] = -Δτ

    # Velocity Constraints
    for i = 1:N
        ψ_vec[1+i] = norm(Δv_vec[i, :]) - 1e-2
    end

    return ψ_vec
end

#============================================================
CALCULATE_MISS:

Description: Calculates miss distance between trajectories

Inputs:
    1. Δτ     - Kepler's Universal Variable
    2. Δv_vec - Matrix of size (N, 3) where each row is the velocity vector at a segment of the trajectory
    3. x̄0_P   - Non-dimensional initial state vector
    4. x̄f_E   - Non-dimensional final state vector
    5. N      - Number of segments of the trajectory

Outputs:
    1. Δxf - Miss distance
============================================================#

function calculate_miss(Δτ, Δv_vec, x̄0_P, x̄f_E, N)

    # Propagating To Final State
    xf_P, Δt = prop_stateUV_Nseg(x̄0_P, Δv_vec, Δτ, N)

    # Finding Moon at Final State
    xf_E = propagate_x(x̄f_E, Δt, 1.0)

    # Finding Miss Distance
    Δxf = abs.(xf_E[1:3] - xf_P[1:3])

    return Δxf
end
