#============================================================
CREATE_MOON:

Description: Creates moon object with appropiate vertex and face Values

Inputs:
    1. Name - String of the name of the moon
    2. Radius - Radius of the moon
    3. μ - Gravitational parameter of the moon
    4. FacetValues - Values of each of the moon's Facets
    5. Visits - Faces of the moon that have already been visited

Outputs:
    1. Moon object of interest
============================================================#

function Create_Moon(Name::String, 
    Radius::Float64, 
    μ::Float64,
    FacetValues::Vector{Int},
    Visits::Union{Vector{Int}, Nothing} = nothing)
  
    # Getting Vertex Data
    vertices = read_gtoc_data()
    combinations = read_face_data()
  
    # Useful Constant
    ϕ = (1 + √5)/2      # Golden Ratio
    nϕ = norm([3*ϕ, 1]) # Vector Length
  
    # Generated Sets
    Facets = LazySets.VPolygon[]
    FacetID = Int[]
    for (i, set) in enumerate(combinations)
      if i ∈ [3, 4, 5, 6] # Longitude Edges
  
      elseif i ∈ [18, 19, 24, 25] # Latitude Edges
        points = [[vertices.long[k], vertices.lat[k]] for k in set]
        s = (sign(points[3][1]), sign(points[3][2]))
        push!(points, [0,        s[2]*90])
        push!(points, [s[1]*180, s[2]*90])
        push!(points, [s[1]*180, s[2]*asind(3*ϕ/nϕ)])
        
        push!( Facets, VPolygon(points) )
        push!( FacetID, i )
      else
        points = [[vertices.long[k], vertices.lat[k]] for k in set]
        push!( Facets, VPolygon(points) )
        push!( FacetID, i )
      end
    end
  
    return Moon(Name, Radius, μ, Facets, FacetID, FacetValues, nothing)
end

#============================================================
CREATE_JUPITERMOON:

Description: Function to specifically create a moon object of one of the moons of Jupiter with predefined face Values

Inputs:
    1. Name - Name of moon of interest
    2. C - Planetary Constants Model of moon of interest
    3. Visits - Faces of moon that have been visited

Outputs:
    1. Moon object of interest
============================================================#

function Create_JupiterMoon(Name::String, C::PlanetaryConstantsModel, Visits::Union{Nothing, Vector{Int}}=nothing)
  
    # Checking Name
    sym = Symbol(Name)
    @assert sym ∈ [:Io, :Europa, :Ganymede, :Callisto]
  
    # Creating Face Values
    vals = zeros(Int, 32)
    if sym ∈ [:Io, :Europa]
      vals[1:8] .= 1
      vals[cat(9:14, 27:32, dims=1)] .= 2
      vals[15:26] .= 3
  
    else
      vals[1:8] .= 3
      vals[cat(9:14, 27:32, dims=1)] .= 2
      vals[15:26] .= 1
  
    end
  
    return Create_Moon(Name, C.R, C.μ, vals, Visits)
end