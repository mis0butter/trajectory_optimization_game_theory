module trajectory_optimization_game_theory

using Symbolics: Symbolics, @variables
using ParametricMCPs: ParametricMCPs, ParametricMCP
using BlockArrays: BlockArray, Block, mortar, blocks
using LinearAlgebra: norm_sqr, cross, dot, norm
using LazySets: LazySets
using TrajectoryGamesBase: TrajectoryGamesBase, AbstractDynamics

# # Write your package code here.
include("propagator.jl") 

# include("parametric_game.jl")
# export ParametricGame, total_dim, solve

# include("environment.jl")
# export get_constraints, CubeEnvironment

# end
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

end # module UpdatedAdversarialTourDesign
