#============================================================
READ_GTOC_DATA:

Description: Reads gtoc data from a csv file and converts it to Latitude and Longitude coordinates

Inputs:
  1. N/A

Outputs:
  2. Array of vertices in x,y,z and Lat/Long coords
============================================================#

function read_gtoc_data()
    # Pulling Data From CSV
    header = [:x, :y, :z]
    vdata = CSV.File("src/modeling/gtoc6.txt"; delim=" ", header=header)
    # @info "Printing" vdata
    pdata = str -> eval(Meta.parse(replace(str, "p" => "(1 + √5)/2")))
    data = [pdata.(getproperty(vdata, h)) for h in header]
    data = Float64.(hcat(data...))
    
    # Converting to Lat/Long
    Φ = zeros(60)
    λ = zeros(60)
    for i = 1:60
      Φ[i] = atand(data[i, 3]/norm(data[i, :]))
      λ[i] = atand(data[i, 2], data[i, 1])
    end
    
    return (; x = data[:, 1], y = data[:, 2], z = data[:, 3], lat = Φ, long = λ)
end
#============================================================
READ_FACE_DATA:

Description: Reads data from a csv file and creates sets of combinations of faces

Inputs:
  1. N/A

Outputs:
  1. out - sets of combinations of faces
============================================================#
function read_face_data()
  # Pulling Data From CSV
  header = [:first, :second, :third, :forth, :fifth, :sixth]
  vdata = @suppress_err CSV.File("src/modeling/gtoc6_faces.txt"; delim=" ", header=header)
  pdata = str -> str
  data = [pdata.(getproperty(vdata, p)) for p in header]

  # Converting to Sets
  out = Vector{Int}[]
  for i = 1:length(data[1])
    push!(out, filter(!ismissing, [data[j][i] for j in 1:6]))
  end

  return out
end