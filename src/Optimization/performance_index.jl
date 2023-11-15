#============================================================
PERFORMANCE_INDEX:

Description: Creates the performance index from an Augmented Lagrangian and Constraints

Inputs:
    1. x - state vector of the form [Δτ, Δv₁, Δv₂, ..., Δvₙ]
    2. λ - Lagrange Mulitplier
    3. p - Penalty vector
    4. N - Number of segments the trajectory is split into
    5. x̄ₒ - Non-dimensional initial state vector
    6. x̄fₒ - Non-dimensional final state vector
    7. mC - Number of terminal constraints 

Outputs:
    1. ϕ - Performance index
============================================================#

function performance_index(
    x, 
    λ, 
    p, 
    N, 
    x̄₀, 
    x̄f₀, 
    mC
)

    # Unwrapping State
    Δτ = x[1]
    Δv_vec = reshape(collect(x[2:end]), N, 3)

    # Objective
    obj = fuel_norm(Δv_vec, N)

    # Constraints
    h_vec = equality_constraints(Δτ, Δv_vec, x̄₀, x̄f₀, N)
    ψ_vec = inequality_constraints(Δτ, Δv_vec, x̄₀, x̄f₀, N)

    # Constructing Augmented Lagrangian
    ϕ = obj

    # Adding Equality Constraints
    for i = 1:mC
        if iszero(p[i]); continue; end

        dϕ = λ[i]*h_vec[i] + 0.5*p[i]*h_vec[i]^2

        @assert dϕ > 0
        ϕ += dϕ
    end

    # Adding Inequality Constraints
    for i in eachindex(ψ_vec)
        if iszero(p[mC+i]); continue; end

        if ψ_vec[i] > -λ[mC+i]/p[mC+i]

            dϕ = λ[mC+i]*ψ_vec[i] + 0.5*p[mC+i]*ψ_vec[i]^2

            ϕ += dϕ

        else
            ϕ += -λ[mC+i]^2/(2*p[mC+i])

        end
    end

    return ϕ
end

#============================================================
FUEL_NORM:

Description: Helper function for performance_index() that calculates the fuel norm

Inputs:
    1. Δv_vec - Matrix of size (N, 3) where each row is the velocity vector at a segment of the trajectory
    2. N - number of segments of the trajectory

Output:
    1. fuel_norm - the fuel norm
============================================================#

function fuel_norm(Δv_vec, N)

    Δv = [norm(Δv_vec[i, :]) for i in 1:N]

    fuel_norm = sum(Δv)

    return fuel_norm
end

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
    N ) 

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
    N ) 

    ψ_vec = zeros( N+1 )

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
    xf_E = propKepΔt(x̄f_E, Δt, 1.0)

    # Finding Miss Distance
    Δxf = abs.(xf_E[1:3] - xf_P[1:3])

    return Δxf
end
