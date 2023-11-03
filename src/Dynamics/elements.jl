#============================================================

OE2RV:

Description: Converts orbital elements to radius and velocity vectors

Inputs:
    1. oe - Keplerian/AbstractVector containing the orbital elements
    2. μ - Gravitational parameter

Outputs:
    1. r_vec - Radius vector
    2. v_vec - Velocity Vector

============================================================#

function oe2rv(oe::AbstractVector{T}, μ::T = 1.0) where T<:Real
    # Decomposing Vector
    a, e, i, ω, Ω, M = oe

    # Handling Options
    ν = M2ν(M, e; units = π/180)
    cν, sν = cosd(ν), sind(ν)

    # Calculating Integrals of Motion
    r = a*(1-e^2) / (1 + e*cν)
    p = a*(1-e^2);
    v = sqrt(μ / p)

    # Perifocal States
    r_vec′ = r * [cν, sν, 0.0]
    v_vec′ = v * [-sν, (e + cν), 0.0]

    # Rotating to Intertial Coordinates from Perifocal
    R = RotZXZ(Ω, i, ω)

    # Creating Intertial State
    r_vec = R*r_vec′
    v_vec = R*v_vec′

    return r_vec, v_vec
end

#============================================================

RV2OE: 

Description: Converts from radius and velocity vectors to orbital elements

Inputs:
    1. r_vec - Radius vector
    2. v_vec - Velocity vector
    3. μ - Gravitational parameter

Outputs:
    1. oe - Orbital elements

============================================================#

function rv2oe(
    rv_vec,  
    μ = 1.0
    ) 

    r_vec = rv_vec[1:3] 
    v_vec = rv_vec[4:6] 

    # Calculating Norms
    r = norm(r_vec)
    v = norm(v_vec)

    # Accounting for planar positions
    if r_vec[3] == 0
        r_vec[3] = 1e-12
    end

    # Calculating Integrals of Motion
    ẑ = [0., 0., 1.] 

    h_vec = r_vec×v_vec;           
    h = norm(h_vec); 
    ĥ = h_vec/h

    n_vec = ẑ×h_vec            
    n = norm(n_vec)  
    n̂ = n_vec/n

    e_vec = (v_vec×h_vec)/μ - r_vec/r 
    e = norm(e_vec) 
    ê = e_vec/e
    ε = 0.5*v^2 - μ/r

    # Semi-Major Axis
    a = -μ/(2*ε)

    # Inclination
    i = acos(ĥ[3])

    # Longitude of Ascending Node
    Ω = atan(n̂[2], n̂[1])

    # Argument of Peri
    ω = atan((n̂×ê)'*ĥ, (n̂⋅ê))

    # True Anomaly
    ν = atan((ê×r_vec/r)'*ĥ, (ê⋅r_vec/r))

    # Converting to Requested Anomaly
    α = ν2M(ν, e; units = π/180)

    oe = [a, e, i, ω, Ω, α]

    return oe
end 

#============================================================

KEP2CART:

Description: Converts from a state vector to orbital elements

Inputs:
    1. kep - Keplerian object of interest
    2. μ - Gravitational parameter

Outputs:
    1. Cartesian object containing r and v vectors

============================================================#

# function kep2cart(kep::Keplerian{T}, μ::T) where T<:Real
#     r, v = oe2rv(collect(kep), μ)
#     # return Cartesian(r, v, kep.Epoch)
#     return [r, v]
# end

#============================================================

CART2KEP:

Description: Converts from a Cartesian object to Keplerian object

Inputs:
    1. cart - Cartesian object of interest
    2. μ - Gravitational parameter

Outputs:
    1. Keplerian object containing orbital elements

============================================================#

# function cart2kep(cart::Cartesian{T}, μ::T) where T<:Real
#     oe = rv2oe( [cart.Position, cart.Velocity], μ )
#     # return Keplerian(oe..., cart.Epoch)
#     return oe 
# end

#============================================================

PCM2CART:

Description: Converts from a Planetary Constants Model object to a Cartesian object

Inputs:
    1. m - Planetary Constants Model object of interest
    2. μ - Gravitational parameter

Outputs:
    1. Cartesian object containing r and v vectors

============================================================#

# function pcm2cart(m::PlanetaryConstantsModel, μ)
#     r, v = oe2rv(collect(m.Elements),μ)
#     # return Cartesian(r, v, m.Elements.Epoch)
#     return [r, v]
# end

# #============================================================

# X2OE:

# Description: Converts from a state vector to orbital elements

# Inputs:
#     1. x - state vector of the form [r;v]
#     2. μ - Gravitational parameter

# Outputs:
#     1. oe - Orbital elements

# ============================================================#

# function x2oe(x::AbstractVector{T}, μ::T = 1.0) where {T<:Real}
#     r_vec = x[1:3]
#     v_vec = x[4:6]
#     oe = rv2oe(r_vec, v_vec, μ)
#     return oe
# end

# ======================================================= #
# ================= Anomaly Conversions ================= #
# ======================================================= #

#============================================================

M2E:

Description: Converts from Mean Anomoly to Eccentric Anomoly

Inputs:
    1. M - Mean Anomoly
    2. e - Eccentricity
    3. tol - Convergence tolerance, initialized to 1e-6
    4. units - Unit converter, initialized to π/180 to convert from deg to rad

Outputs:
    1. E - Eccentric Anomoly

==============================================================#

function M2E(M, e; tol = 1e-6, units = π/180)
    # Creating initial guess for E
    M *= units
    E = M < π ? M + e/2 : M - e/2

    # Iterating
    ΔE = Inf
    iter = 0
    while abs(ΔE) > tol
        iter += 1
        ΔE = (E - e*sin(E) - M)/(1 - e*cos(E))
        E -= ΔE
    end
    
    return E/units
end

#============================================================

E2M:

Description: Converts from Eccentric Anomoly to Mean Anomoly

Inputs:
    1. E - Eccentric Anomoly
    2. e - Eccentricity
    3. units - Unit converter, initialized to π/180 to convert from deg to rad

Outputs:
    1. M - Mean Anomaly

============================================================#

function E2M(E, e; units = π/180)
    M = (E*units - e*sin(E*units))/units
    return M
end

#============================================================

E2ν:

Description: Converts from Eccentric Anomaly to True Anomaly

Inputs:
    1. E - Eccentric Anomoly
    2. e - Eccentricity
    3. units - Unit converter, initialized to π/180 to convert from deg to rad

Outputs:
    1. ν - True Anomoly

============================================================#

function E2ν(E, e; units = π/180)
    ν = 2*atan(sqrt((1+e)/(1-e))*tan(E/2*units)) / units 
    return ν
end

#============================================================

ν2E:

Description: Converts from True Anomaly to Eccentric Anomaly

Inputs:
    1. ν - True Anomoly
    2. e - Eccentricity
    3. units - Unit converter, initialized to π/180 to convert from deg to rad

Outputs:
    1. E - Eccentric Anomoly

============================================================#

function ν2E(ν, e; units = π/180)
    E = 2*atan(sqrt((1-e)/(1+e))*tan(ν/2*units)) / units 
    return E
end

#============================================================

ν2M:

Description: Converts from True Anomaly to Mean Anomaly

Inputs:
    1. ν - True Anomoly
    2. e - Eccentricity
    3. units - Unit converter, initialized to π/180 to convert from deg to rad

Outputs:
    1. M - Mean Anomoly

============================================================#

function ν2M(ν, e; units = π/180)
    M = E2M(ν2E(ν, e; units = units), e; units = units)
    return M
end

#============================================================

M2ν:

Description: Converts from Mean Anomaly to True Anomaly

Inputs:
    1. M - Mean Anomoly
    2. e - Eccentricity
    3. units - Unit converter, initialized to π/180 to convert from deg to rad

Outputs:
    1. ν - True Anomoly

============================================================#

function M2ν(M, e; units = π/180)
    ν = E2ν(M2E(M, e; units = units), e; units = units)
    return ν
end 

