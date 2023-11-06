# =====================================================================
# === Lambert

# Main structure
struct Lambert{T}
  xᵢ::AbstractVector{T}
  xⱼ::AbstractVector{T}
  tof::T
  Δvᵢ::AbstractVector{T}
  Δvⱼ::AbstractVector{T}
  _iterations
  _error
  _code
end

include("lambert.jl")

# Exports
export lambert

# =====================================================================
# === Moon

# Main structure
struct Moon
    Name::String
    Radius::Float64
    μ::Float64
    Facets::Vector{LazySets.VPolygon}
    FacetID::Vector{Int}
    FacetValues::Vector{Int}
    Visits::Union{Vector{Int}, Nothing}
  end

include("moon.jl")

# Exports
export Create_JupiterMoon, Create_Moon

# =====================================================================
# === Project Constants
include("constants.jl")

# Exports
export PlanetaryConstantsModel, ProblemConstants, import_constants

# =====================================================================
# === Plotting
include("plotting.jl")

# Exports
export plot_solution!, plot_moon!, plottraj

# =====================================================================
# === Parsing Data
include("parse_data.jl")

# Exports
export read_gtoc_data, read_face_data


