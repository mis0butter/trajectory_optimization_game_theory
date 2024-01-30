module trajectory_optimization_game_theory

# Math
using LinearAlgebra: norm 
using StaticArrays          # don't think I use this 
using Rotations: AngleAxis, RotZ, RotZXZ

# Optimization
using Roots: find_zero      # don't think I use this  
using ForwardDiff           # think I only use gradient 
using Optim                 # think i only use ... some functions 

# Modeling
using LazySets              # don't think I use this 
using JuMP: JuMP, @variable, @constraint, @objective
using OSQP: OSQP 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

# Tools
using CSV
using Suppressor: @suppress_err

## ============================================ ##

# Junette and Sofia Dyn 
include("Dyn/Dyn.jl") 
include("Opt/Opt.jl") 
include("Utils/Utils.jl")
include("Lambert/Lambert.jl") 
include("Opt/matrix_game_solver.jl")

## ============================================ ##
#  Sub-Module Includes

# include("old/Dynamics/Dynamics.jl")
# include("old/Modeling/Modeling.jl")
# include("old/Optimization/Optimization.jl")


end # module UpdatedAdversarialTourDesign


