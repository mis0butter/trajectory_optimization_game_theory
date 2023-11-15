#============================================================
NONDIMENSIONALIZE_X:

Description: Create relative distance and time constants and non-dimensionalize a state vector.

Inputs:
    1. x_vec - state vector of the form [r;v]
    2. m - Planetary Constants Model object 
    3. R - Radius of the planet 

Outputs:
    1. x̄_vec - Non-dimensionalized x_vec
    2. DU - Distance unit
    3. TU - Time unit
============================================================#

function nondimensionalize_x(
    x_vec, 
    μ,
    R )

    # Separating State Vectors
    r_vec = x_vec[1:3]
    v_vec = x_vec[4:6]

    # Passing to Base Routine
    r_vec, v_vec, DU, TU = nondimensionalize_rv(r_vec, v_vec, μ, R)
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
    4. R - radius of the planet 

Outputs:
    1. r̄_vec - Non-dimensionalized radius vector
    2. v̄_vec - Non-dimensionalized velocity vector
    3. DU - Distance unit
    4. TU - Time unit
============================================================#

function nondimensionalize_rv(
    r_vec, 
    v_vec, 
    μ,
    R ) 

    # Distance unit DU is defined by the Earth radius 
    DU = R

    # Time unit TU is defined by Earth mu
    TU = sqrt(DU^3/μ)

    # Converting Units
    r̄_vec = r_vec / DU
    v̄_vec = v_vec / (DU/TU)

    # Outputting
    return r̄_vec, v̄_vec, DU, TU
end