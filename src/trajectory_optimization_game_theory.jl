using Infiltrator 
using ForwardDiff 
using LinearAlgebra 

module trajectory_optimization_game_theory

# Write your package code here.
include("Dynamics/Dynamics.jl") 

include("Optimization/Optimization.jl") 

function test_internal(x) 
    return 2*x 
end 

function test_external(x) 
    return 2*test_internal(x) 
end 
export test_external 

end
