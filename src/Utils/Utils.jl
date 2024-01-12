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

## ============================================ ##

""" 
Get orthogonal axes of local frame at given state vector: 
    axis_1: along velocity vector 
    axis_2: along radius vector 
    axis_3: normal to orbit plane 
""" 
function axis_123( rv_vec ) 

    # center of polygon 
    r_f = rv_vec[1:3]  ; axis_2 = -r_f / norm(r_f)
    v_f = rv_vec[4:6]  ; axis_1 = v_f / norm(v_f)   

    # define vector normal to orbit plane 
    axis_3 = cross( axis_1, axis_2 )  ; axis_3 = axis_3 / norm(axis_3) 

    return axis_1, axis_2, axis_3 
end 

export axis_123 

## ============================================ ##

"Compute vertices of polygon around given state vector"  
function polygon_vertices( 
    rv_vec,                 # [N,6] state vector 
    dist = 6378.0 / 10,    # radius of polygon 
) 

    # center of polygon 
    r_vec = rv_vec[1:3] 
    _, axis_2, axis_3 = axis_123( rv_vec ) 

    # top vertex: move up from r_f along axis 3 
    r_top = r_vec + axis_3 * dist 

    # top-inner vertex: move up from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = cosd(60) * axis_3 * dist + sind(60) * axis_2 * dist 
    r_topin  = r_vec + vec

    # bottom-inner vertex: move down from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * dist + sind(60) * axis_2 * dist 
    r_botin  = r_vec + vec 

    # bottom vertex: move down from r_f along axis 3 
    r_bot = r_vec - axis_3 * dist 

    # bottom-outer vertex: move down from r_f along axis 3 and right along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * dist - sind(60) * axis_2 * dist 
    r_botout = r_vec + vec 

    # top-outer vertex: move up from r_f along axis 3 and right along axis 2, 60 degrees 
    vec       = cosd(60) * axis_3 * dist - sind(60) * axis_2 * dist 
    r_topout  = r_vec + vec 

    vertices = ( top = r_top, topin = r_topin, botin = r_botin, bot = r_bot, botout = r_botout, topout = r_topout ) 
    return vertices 
end 

