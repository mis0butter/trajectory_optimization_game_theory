# =====================================================================
# === Broyden–Fletcher–Goldfarb–Shanno

# Including Methods
include("bfgs.jl")

# Exports
export BFGS, bfgs

# =====================================================================
# === Sims-Flanagan

# Main Structure
struct SimsFlanaganTrajectory{T<:AbstractFloat}
    Δτ::T
    Δv⃗::AbstractMatrix{T}
    Δv̄::AbstractMatrix{T}
    Δv::T
    DU::T
    TU::T
    x0::AbstractVector{T}
    xf::AbstractVector{T} 
    Δt::T
    miss::AbstractVector{T}
end

# Including Methods
include("performance_index.jl")
include("constraints.jl")
include("nondimm.jl")
include("state_update.jl")
include("solve_transfer.jl")

# Exports
# export SimsFlanaganTrajectory, solve_transfer