"""
    Replaces environment.jl in TrajectoryGamesBase for a SphereEnvironment instead of a PolygonEnvironment
"""

function get_constraints end
abstract type AbstractEnvironment end

struct CubeEnvironment{T} <: AbstractEnvironment
    set::T
end

function CubeEnvironment(sideLength::Float64 = 1.0)
     
    # define unit cube
    x = [0.5, 0.5, 0.5, 0.5, -0.5, -0.5, -0.5, -0.5]
    y = [-0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5, -0.5]
    z = [-0.5, -0.5, 0.5, 0.5, -0.5, -0.5, 0.5, 0.5]

    # stack coordinates and scale by input dimension sideLength
    bounds = [x y z] * sideLength

    # convert matrix into vector of vectors 
    bounds = [bounds[i, :] for i in 1:size(bounds,1)]
    CubeEnvironment(bounds)
end

function CubeEnvironment(vertices::AbstractVector{<:AbstractVector{<:Real}})
    CubeEnvironment(LazySets.VPolytope(vertices))
end

function get_constraints(environment::CubeEnvironment, player_index = Nothing)
    constraints = LazySets.constraints_list(environment.set)
    function (state)
        positions = (substate[1:3] for substate in blocks(state))
        mapreduce(vcat, Iterators.product(constraints, positions)) do (constraint, position)
            -constraint.a' * position + constraint.b
        end
    end
end