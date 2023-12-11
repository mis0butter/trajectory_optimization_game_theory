include("plotting.jl")
include("IC.jl") 

## ============================================ ##

"Convert vector of vectors into matrix"
function vv2m( vec_vec ) 

    out = mapreduce( permutedims, vcat, vec_vec ) 

    return out 
end 

export vv2m 

## ============================================ ##

"Return index that matches val"
function get_index( 
    A,      # array (or vector)  
    val,    # value to find 
) 

    ind = findall( x -> x == val, A )

    return ind 
end 

export get_index 

## ============================================ ##

"Convert matrix into vector of vectors" 
function m2vv( M )

    N       = size(M, 1) 
    vec_vec = [] 
    for i = 1 : N 
        push!( vec_vec, M[i,:] ) 
    end 

    return vec_vec 
end 

export m2vv 