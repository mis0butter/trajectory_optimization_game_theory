include("plotting.jl")
include("IC.jl") 

## ============================================ ##

"Convert vector of vectors into matrix"
function vv2m( vec_vec ) 

    out = mapreduce( permutedims, vcat, vec_vec ) 

    return out 
end 

export vv2m 

