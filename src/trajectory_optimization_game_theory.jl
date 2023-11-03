module trajectory_optimization_game_theory

using Symbolics: Symbolics, @variables
using ParametricMCPs: ParametricMCPs, ParametricMCP
using BlockArrays: BlockArray, Block, mortar, blocks
using LinearAlgebra: norm_sqr
using LazySets: LazySets

# Write your package code here.
include("propagator.jl") 

include("parametric_game.jl")
export ParametricGame, total_dim, solve

include("environment.jl")
export get_constraints, CubeEnvironment

end
