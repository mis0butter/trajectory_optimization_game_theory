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
# === Plotting
include("plotting.jl")

# Exports
export plot_solution!, plot_sims_flanagan!, plottraj



