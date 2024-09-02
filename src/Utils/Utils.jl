include("plotting.jl")
include("IC.jl") 
include("structs.jl")

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
    dist = R_polygon,    # radius of polygon 
) 

    # center of polygon 
    r_vec    = rv_vec[1:3] 
    _, axis_2, axis_3 = axis_123( rv_vec ) 

    # top vertex: move up from r_f along axis 3 
    r_top    = r_vec + axis_3 * dist 

    # top-inner vertex: move up from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = cosd(60) * axis_3 * dist + sind(60) * axis_2 * dist 
    r_topin  = r_vec + vec

    # bottom-inner vertex: move down from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * dist + sind(60) * axis_2 * dist 
    r_botin  = r_vec + vec 

    # bottom vertex: move down from r_f along axis 3 
    r_bot    = r_vec - axis_3 * dist 

    # bottom-outer vertex: move down from r_f along axis 3 and right along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * dist - sind(60) * axis_2 * dist 
    r_botout = r_vec + vec 

    # top-outer vertex: move up from r_f along axis 3 and right along axis 2, 60 degrees 
    vec      = cosd(60) * axis_3 * dist - sind(60) * axis_2 * dist 
    r_topout = r_vec + vec 

    vertices = ( top = r_top, topin = r_topin, botin = r_botin, bot = r_bot, botout = r_botout, topout   = r_topout ) 
    return vertices 
end 

export polygon_vertices 

## ============================================ ##

"Generate uniformly distributed random point(s) within a circle of radius R"

function unif_random_points_circle( 
    R = R_polygon,  # radius of circle 
    N = 1,          # number of points 
) 

    θ = 2*pi*rand(N) 
    r = R * sqrt.(rand(N)) 

    x = r .* cos.(θ) 
    y = r .* sin.(θ) 

    return x, y 
end 

export unif_random_points_circle 

## ============================================ ##

"Generate uniformly distributed random point(s) within a circle of radius R around given state vector" 

function rand_IC( 
    rv_vec,         # [N,6] state vector 
    R = R_polygon,  # radius of circle 
    N = 1,          # number of points 
) 

    # center of polygon 
    r_vec = rv_vec[1:3] 
    _, axis_2, axis_3 = axis_123( rv_vec ) 

    # generate random points 
    x, y = unif_random_points_circle( R, N ) 

    vec   = x[1] .* axis_2 + y[1] .* axis_3 
    r_out = r_vec + vec 

    return r_out 
end

export rand_IC 

## ============================================ ##

"Compute norm of U vectors for both players"
function U_norm( game, params ) 

    k_tt_replan = params.k_tt_replan 
    k_max       = game.k_replan[end] 

    p1_U_hist = [] 
    p2_U_hist = [] 
    for kk = 1 : k_max 

        p1_chosen, p2_chosen = p_strategy( game, kk, params.strategy ) 

        p1_U = game.p1_state[ kk ].U[ p1_chosen ][ 1 : k_tt_replan, : ] 
        p2_U = game.p2_state[ kk ].U[ p2_chosen ][ 1 : k_tt_replan, : ] 
    
        push!( p1_U_hist, p1_U ) 
        push!( p2_U_hist, p2_U ) 
    end 
    
    p1_U_hist = mapreduce( permutedims, hcat, p1_U_hist )' 
    p2_U_hist = mapreduce( permutedims, hcat, p2_U_hist )' 
    
    p1_U_norm = [ norm( p1_U_hist[ii,:] ) for ii in 1 : size(p1_U_hist, 1) ] 
    p2_U_norm = [ norm( p2_U_hist[ii,:] ) for ii in 1 : size(p2_U_hist, 1) ] 
    
    return p1_U_norm, p2_U_norm 
end 

export U_norm 


## ============================================ ##

function p_rv_ref_hist( game, params ) 

    k_tt_replan = params.k_tt_replan 
    k_max       = game.k_replan[end] 

    tt_hist     = [] 
    p1_rv_hist  = [] 
    p2_rv_hist  = [] 
    rv_ref_hist = [] 
    for kk = 1 : k_max 

        p1_chosen, p2_chosen = p_strategy( game, kk, params.strategy ) 
    
        # get traveled trajectory 
        tt  = game.t_ref_E[ kk ][ 1 : k_tt_replan ] 
        p1_r = game.p1_state[ kk ].X[ p1_chosen ][ 1 : k_tt_replan , : ] 
        p2_r = game.p2_state[ kk ].X[ p2_chosen ][ 1 : k_tt_replan , : ] 
        rv_ref = game.rv_ref_E[ kk ][ 1 : k_tt_replan , : ] 
        if kk == k_max 
            tt   = game.t_ref_E[ kk ][ 1 : k_tt_replan + 1 ] 
            p1_r = game.p1_state[ kk ].X[ p1_chosen ][ 1 : k_tt_replan + 1, : ] 
            p2_r = game.p2_state[ kk ].X[ p2_chosen ][ 1 : k_tt_replan + 1, : ] 
            rv_ref = game.rv_ref_E[ kk ][ 1 : k_tt_replan + 1, : ] 
        end 

        push!( tt_hist, tt ) 
        push!( p1_rv_hist, p1_r ) 
        push!( p2_rv_hist, p2_r ) 
        push!( rv_ref_hist, rv_ref ) 

    end 
    
    tt_hist     = mapreduce( permutedims, hcat, tt_hist )' 
    p1_rv_hist  = mapreduce( permutedims, hcat, p1_rv_hist )' 
    p2_rv_hist  = mapreduce( permutedims, hcat, p2_rv_hist )' 
    rv_ref_hist = mapreduce( permutedims, hcat, rv_ref_hist )' 

    @exfiltrate 

    return tt_hist, p1_rv_hist, p2_rv_hist, rv_ref_hist  
end 

export p_rv_ref_hist 


## ============================================ ## 

function dist_ref_norm( game, params ) 

    _, p1_rv_hist, p2_rv_hist, rv_ref_hist = p_rv_ref_hist( game, params )

    p1_r_diff = p1_rv_hist[:,1:3] - rv_ref_hist[:,1:3] 
    p2_r_diff = p2_rv_hist[:,1:3] - rv_ref_hist[:,1:3] 

    p1_ref_norm = [ norm( p1_r_diff[ii,:] ) for ii in 1 : size(p1_r_diff, 1) ] 
    p2_ref_norm = [ norm( p2_r_diff[ii,:] ) for ii in 1 : size(p2_r_diff, 1) ] 
    
    return p1_ref_norm, p2_ref_norm 
end 

export dist_ref_norm 


## ============================================ ## 

function dist_norm( game, params ) 

    _, p1_rv_hist, p2_rv_hist, _ = p_rv_ref_hist( game, params )

    r_diff = p1_rv_hist[:,1:3] - p2_rv_hist[:,1:3] 
    r_norm = [ norm( r_diff[ii,:] ) for ii in 1 : size(r_diff, 1) ] 
    
    return r_norm
end 

export dist_norm 

