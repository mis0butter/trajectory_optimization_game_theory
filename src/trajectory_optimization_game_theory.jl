module trajectory_optimization_game_theory

# Write your package code here.
include("propagator.jl") 

## ============================================ ##

# module AdversarialSpacecraftTrajectories

# Math
using LinearAlgebra
using StaticArrays
using Rotations: AngleAxis, RotZ, RotZXZ

# Optimization
using Roots: find_zero
using ForwardDiff

# Modeling
using LazySets

# Tools
using CSV
using Suppressor: @suppress_err

# =====================================================================
# === Sub-Module Includes

include("Dynamics/Dynamics.jl")

include("Modeling/Modeling.jl")

include("Optimization/Optimization.jl")

## ============================================ ##

include("Lambert/Lambert.jl") 

end # module UpdatedAdversarialTourDesign


