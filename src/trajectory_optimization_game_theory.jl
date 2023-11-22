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
include("Opt/Opt.jl") 
include("Utils/Utils.jl")
include("Lambert/Lambert.jl") 

## ============================================ ##
#  Sub-Module Includes

include("old/Dynamics/Dynamics.jl")
include("old/Modeling/Modeling.jl")
include("old/Optimization/Optimization.jl")


end # module UpdatedAdversarialTourDesign


