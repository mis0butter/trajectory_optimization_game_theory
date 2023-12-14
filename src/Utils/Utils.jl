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

    ind = findall( x -> x == val, A )[1] 

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

## ============================================ ##

"Compute Δv from desired inclination change" 
function computeInclinationChange(
    rv,     # initial state vector 
    Δi,     # desired inclination change 
    mu,     # gravitational parameter 
)

    v    = rv[4:6]
    kep2 = cart2kep(rv, mu)

    kep2[3] += Δi
    v_des    = kep2cart(kep2, mu)[4:6]

    # compute Δv from geometry 
    Δi_mag = sqrt(norm(v)^2 + norm(v_des)^2 - 2*norm(v)*norm(v_des)*cos(Δi))
    Δv = Δi_mag * (v_des - v) 

    return Δv 
end

export computeInclinationChange

