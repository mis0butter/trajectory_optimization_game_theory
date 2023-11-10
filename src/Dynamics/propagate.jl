#============================================================
PROPAGATE_PLANETARYCONSTANTSMODEL:

Description: Solves for Mean Anomoly/Motion of a Planetary Constants Model

Inputs:
    1. m - Planetary Constants Model of interest
    2. Δt - Time period
    3. μ - Gravitational parameter

Outputs:
    1. Planetary Constants Model of interest at end of time period
============================================================#

function propagate_PlanetaryConstantsModel(m::PlanetaryConstantsModel{T}, Δt::Real, μ::T)::PlanetaryConstantsModel{T} where {T<:Real}
    return PlanetaryConstantsModel(m.μ, m.R, m.n, propagate_Keplerian(m.Elements, Δt, μ))
end

#============================================================
PROPAGATE_KEPLERIAN:

Description: Solves for Mean Anomoly/Motion of a Keplerian object at the end of the time period.

Inputs:
    1. K - Keplerian object of interest
    2. Δt - Time period
    3. μ - Gravitational parameter

Outputs:
    1. Keplerian object at the end of the time period
============================================================#

function propagate_Keplerian(K::Keplerian{T}, Δt::Real, μ::T)::Keplerian{T} where {T<:Real}
    # Unwrapping Elements
    oe = [K.SemiMajorAxis, K.Eccentricity, K.Inclination, 
        K.LongitudeAscNode, K.ArgOfPeriapsis, K.MeanAnomaly]

    # Calculating Mean Motion in Deg/sec
    n = rad2deg(sqrt(μ/oe[1]^3))

    # Updating
    oe[6] += n*Δt

    return Keplerian(oe..., K.Epoch + Δt) 
end

#============================================================
PROPAGATE_X:

Description: Solves for the final state vector at the end of a time period

Inputs: 
    1. x₀_vec - Initial state vector
    2. Δt - Time period
    3. μ - Gravitational parameter

Outputs:
    1. Y - Final state vector
============================================================#

function propagate_x(x₀_vec::AbstractVector{T2}, Δt::T1, μ::T2) where {T1<:Real, T2<:Real}

    # len = length(Δt)
    Y = zeros(T1, 6)
    # t = zeros(T1)

    a = -0.5*μ/(0.5*norm(x₀_vec[4:6])^2 - μ/norm(x₀_vec[1:3]))
    # println(a)
    if a > 0
        Y = propKepTE(x₀_vec, Δt, μ)
    else
        Y = propKepTH(x₀_vec, Δt, μ)
    end

    return Y
end

#============================================================
PROPKEPTE:

Description: Elliptical propagation algorithm using Kepler's equations

Inputs:
    1. x₀_vec - Initial state vector
    2. Δt - Time period
    3. μ - Gravitational parameter

Outputs:
    1. xf_vec - Final state vector
============================================================#

function propKepTE(
    x₀_vec::AbstractVector, 
    Δt,
    μ
    )

    # Defining Constants
    r_vec = x₀_vec[1:3]; r = norm(r_vec) 
    v_vec = x₀_vec[4:6]; v = norm(v_vec)
    ε = (0.5*v^2 - μ/r)
    a = -0.5*μ/ε 
    σ = dot(r_vec,v_vec)/sqrt(μ) # Frequent Constant

    # Root Finding for Change in Eccentric Anomaly
    ΔEᵢ = Δt*sqrt(μ/a^3)
    func = ΔE -> -ΔEᵢ + ΔE + (σ/sqrt(a))*(1 - cos(ΔE)) - (1 - r/a)*sin(ΔE)
    ΔE = find_zero(func, ΔEᵢ)

    # Creating Another Useful Constant
    ρ = a + (r - a)*cos(ΔE) + σ*(sqrt(a))*sin(ΔE)

    # Creating f, g, ḟ, ġ
    f = 1 - (a / r) * (1 - cos(ΔE))
    g = a * σ / sqrt(μ) * (1 - cos(ΔE)) + r * sqrt(a / μ) * sin(ΔE)
    ḟ = -sqrt(μ * a) / (ρ * r) * sin(ΔE)
    ġ = 1 - a / ρ * (1 - cos(ΔE))

    # Converting to Cartesian States
    rₖ_vec = f*r_vec + g*v_vec
    vₖ_vec = ḟ*r_vec + ġ*v_vec
    
    xf_vec = vcat(rₖ_vec, vₖ_vec)
    
    return xf_vec
end

#============================================================
PROPKEPTH:

Description: Hyperbolic propagation algorithm using Kepler's equations

Inputs:
    1. x₀_vec - Initial state vector
    2. Δt - Time period
    3. μ - Gravitational parameter

Outputs:
    1. xf_vec - Final state vector
============================================================#

function propKepTH(
    x₀_vec::AbstractVector, 
    Δt,
    μ
    )

    # Defining Constants
    r_vec = x₀_vec[1:3]; r = norm(r_vec) 
    v_vec = x₀_vec[4:6]; v = norm(v_vec)
    ε = (0.5*v^2 - μ/r)
    a = -0.5*μ/ε 
    σ = dot(r_vec,v_vec)/sqrt(μ) # Frequent Constant

    # Root Finding for Change in Hyperbolic Anomaly
    ΔHᵢ = sign(Δt)
    func  = ΔH -> -Δt*sqrt(-μ/a^3) - ΔH + (σ/sqrt(-a))*(cosh(ΔH) - 1) + (1 - r/a)*sinh(ΔH)
    ΔH = find_zero(func, ΔHᵢ)

    # Creating Another Useful Constant
    ρ = a + (r - a)*cosh(ΔH) + σ*(sqrt(-a))*sinh(ΔH)

    # Creating f, g, ḟ, ġ
    f = 1 - (a / r) * (1 - cosh(ΔH))
    g = a * σ / sqrt(μ) * (1 - cosh(ΔH)) + r * sqrt(-a / μ) * sinh(ΔH)
    ḟ = -sqrt(-μ * a) / (ρ * r) * sinh(ΔH)
    ġ = 1 - a / ρ * (1 - cosh(ΔH))

    # Converting to Cartesian States
    rₖ_vec = f*r_vec + g*v_vec
    vₖ_vec = ḟ*r_vec + ġ*v_vec
    
    xf_vec = vcat(rₖ_vec, vₖ_vec)

    return xf_vec
end