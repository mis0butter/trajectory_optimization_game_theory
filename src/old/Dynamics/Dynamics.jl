# =====================================================================
# === Dynamics Basis

abstract type EphemerisState end

struct Cartesian{T} <: EphemerisState 
    # State
    Position::AbstractVector{T}
    Velocity::AbstractVector{T}

    # Time
    Epoch::T
end

struct Keplerian{T} <: EphemerisState 
    # State
    SemiMajorAxis::T
    Eccentricity::T
    Inclination::T
    ArgOfPeriapsis::T
    LongitudeAscNode::T
    MeanAnomaly::T

    # Time
    Epoch::T
end

struct PlanetaryConstantsModel{T}
    # Physical Properties
    μ::T
    R::T
    n::T
  
    # Initial State
    Elements::Keplerian{T}
  end

export EphemerisState, Cartesian, Keplerian, PlanetaryConstantsModel

# Modifying Base Functions
import Base: collect

collect(kep::Keplerian) = [
    kep.SemiMajorAxis,
    kep.Eccentricity,
    kep.Inclination,
    kep.ArgOfPeriapsis,
    kep.LongitudeAscNode,
    kep.MeanAnomaly
]

collect(cart::Cartesian) = vcat(cart.Position, cart.Velocity)


# =====================================================================
# === Keplerian Elements
include("elements.jl")

# Exports
export rv2oe, oe2rv, pcm2cart, cart2kep, kep2cart


# =====================================================================
# === Time Propagation
include("propagate.jl")

export propagate_PlanetaryConstantsModel, propagate_Keplerian, propagate_x

# =====================================================================
# === Universal Variable Propagation
include("universal.jl")

# Exports
export propKepUV