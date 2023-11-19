module trajectory_optimization_game_theory

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

## ============================================ ##

# Junette and Sofia Dyn 
include("Dyn/Dyn.jl") 
include("Utils/Utils.jl")

# =====================================================================
# === Sub-Module Includes

include("Dynamics/Dynamics.jl")

include("Modeling/Modeling.jl")

include("Optimization/Optimization.jl")

## ============================================ ##

include("Lambert/Lambert.jl") 

end # module UpdatedAdversarialTourDesign


