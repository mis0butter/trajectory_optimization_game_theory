#============================================================
NONDIMENSIONALIZE_X:

Description: Create relative distance and time constants and non-dimensionalize a state vector.

Inputs:
    1. x_vec - state vector of the form [r;v]
    2. μ - gravitational parameter

Outputs:
    1. x̄_vec - Non-dimensionalized x_vec
    2. DU - Distance unit
    3. TU - Time unit
============================================================#

function nondimensionalize_x(x_vec::AbstractVector{<:Real}, μ::Real)
    # Separating State Vectors
    r_vec = x_vec[1:3]
    v_vec = x_vec[4:6]

    # Passing to Base Routine
    r_vec, v_vec, DU, TU = nondimensionalize_rv(r_vec, v_vec, μ)
    x̄_vec = vcat(r_vec, v_vec)
    
    # Passing Results
    return x̄_vec, DU, TU
end

#============================================================
NONDIMENSIONALIZE_RV:

Description: Create distance and time units and non-dimensionalize the radius and velocity vectors

Inputs:
    1. r_vec - Radius vector
    2. v_vec - Velocity vectors
    3. μ - Gravitational parameter

Outputs:
    1. r̄_vec - Non-dimensionalized radius vector
    2. v̄_vec - Non-dimensionalized velocity vector
    3. DU - Distance unit
    4. TU - Time unit
============================================================#

function nondimensionalize_rv(r_vec::V, v_vec::V, μ::Real) where V<:AbstractVector{<:Real}
    # Distance unit
    DU = norm(r_vec)

    # Time unit
    vc = sqrt(μ/DU)        # Circular Velocity
    TU = DU/vc

    # Converting Units
    r̄_vec = r_vec / DU
    v̄_vec = v_vec / (DU/TU)

    # Outputting
    return r̄_vec, v̄_vec, DU, TU
end

