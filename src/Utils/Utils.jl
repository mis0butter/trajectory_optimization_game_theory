include("plotting.jl")
include("IC.jl") 
include("structs.jl")

## ====================================================================

"Convert vector of vectors into matrix"
function vv2m( vec_vec ) 

    out = mapreduce( permutedims, vcat, vec_vec ) 

    return out 
end 

export vv2m 

## ====================================================================

"Return index that matches val"
function get_index( 
    A,      # array (or vector)  
    val,    # value to find 
) 

    ind = findall( x -> x == val, A )[1] 

    return ind 
end 

export get_index 

## ====================================================================

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

## ====================================================================

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

## ====================================================================

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

## ====================================================================

"Compute vertices of polygon around given state vector"  
function polygon_vertices( 
    rv_state,                 # [N,6] state vector 
    parameters,    # radius of polygon 
) 

    dist = parameters.R_polygon 

    # center of polygon 
    r_vec    = rv_state[1:3] 
    _, axis_2, axis_3 = axis_123( rv_state ) 

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

## ====================================================================

"Generate uniformly distributed random point(s) within a circle of radius R"

function unif_random_points_circle( 
    R,      # parameters struct 
    rng,    # random number generator 
    N = 1,  # number of points 
) 

    θ = 2*pi*rand(rng, N) 
    r = R * sqrt.(rand(rng, N)) 

    x = r .* cos.(θ) 
    y = r .* sin.(θ) 

    return x, y 
end 

export unif_random_points_circle 

## ====================================================================

"Generate uniformly distributed random point(s) within a circle of radius R around given state vector" 

function rand_IC( 
    rv_reference,   # [N,6] state vector 
    R_polygon,      # parameters struct 
    rng,            # random number generator 
    N = 1,          # number of points 
) 

    # reference vector and frame  
    r_reference = rv_reference[1:3] 
    _, axis_2, axis_3 = axis_123( rv_reference ) 

    # generate random points 
    x_random, y_random = unif_random_points_circle( R_polygon, rng, N ) 

    random_vector = x_random[1] .* axis_2 + y_random[1] .* axis_3 
    r_IC          = r_reference + random_vector 

    return r_IC 
end 

export rand_IC 

## ====================================================================

"Compute norm of U vectors for both players"
function p1_p2_u_hist( game, params = game.params[1] ) 

    k_tt_replan = params.k_tt_replan 
    k_max       = game.k_replan[end] 

    p1_U_hist = [] 
    p2_U_hist = [] 
    for kk = 1 : k_max 

        # propagate SC state forward 
        p1 = game.p1_state[ kk ] 
        p2 = game.p2_state[ kk ] 
        p1_chosen = p1.chosen 
        p2_chosen = p2.chosen 

        p1_U = game.p1_state[ kk ].U[ p1_chosen ][ 1 : k_tt_replan, : ] 
        p2_U = game.p2_state[ kk ].U[ p2_chosen ][ 1 : k_tt_replan, : ] 
    
        push!( p1_U_hist, p1_U ) 
        push!( p2_U_hist, p2_U ) 
    end 
    
    p1_U_hist = mapreduce( permutedims, hcat, p1_U_hist )' 
    p2_U_hist = mapreduce( permutedims, hcat, p2_U_hist )' 
    
    return p1_U_hist, p2_U_hist 
end 

export p1_p2_u_hist 


## ====================================================================

"Compute norm of U vectors for both players"
function U_norm( game, params = game.params[1] ) 

    k_tt_replan = params.k_tt_replan 
    k_max       = game.k_replan[end] 

    p1_U_hist = [] 
    p2_U_hist = [] 
    for kk = 1 : k_max 

        # propagate SC state forward 
        p1 = game.p1_state[ kk ] 
        p2 = game.p2_state[ kk ] 
        p1_chosen = p1.chosen 
        p2_chosen = p2.chosen 

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


## ====================================================================

# tt_hist, p1_rv_hist, p2_rv_hist, rv_ref_hist = p_rv_ref_hist( game, params ) 
function p_rv_ref_hist( game, params = game.params[1] ) 

    k_tt_replan = params.k_tt_replan 
    k_max       = game.k_replan[end] 

    tt_hist     = [] 
    p1_rv_hist  = [] 
    p2_rv_hist  = [] 
    rv_ref_hist = [] 
    for kk = 1 : k_max 

        # propagate SC state forward 
        p1 = game.p1_state[ kk ] 
        p2 = game.p2_state[ kk ] 
        p1_chosen = p1.chosen 
        p2_chosen = p2.chosen 
    
        # get traveled trajectory 
        tt     = game.t_ref_E[ kk ][ 1 : k_tt_replan ] 
        p1_r   = game.p1_state[ kk ].X[ p1_chosen ][ 1 : k_tt_replan , : ] 
        p2_r   = game.p2_state[ kk ].X[ p2_chosen ][ 1 : k_tt_replan , : ] 
        rv_ref = game.rv_ref_E[ kk ][ 1 : k_tt_replan , : ] 
        if kk == k_max 
            tt     = game.t_ref_E[ kk ][ 1 : k_tt_replan + 1 ] 
            p1_r   = game.p1_state[ kk ].X[ p1_chosen ][ 1 : k_tt_replan + 1, : ] 
            p2_r   = game.p2_state[ kk ].X[ p2_chosen ][ 1 : k_tt_replan + 1, : ] 
            rv_ref = game.rv_ref_E[ kk ][ 1 : k_tt_replan + 1, : ] 
        end 

        push!( tt_hist, tt ) 
        push!( p1_rv_hist, p1_r ) 
        push!( p2_rv_hist, p2_r ) 
        push!( rv_ref_hist, rv_ref ) 

    end 
    
    tt_hist     = float.( mapreduce( permutedims, hcat, tt_hist )'     )[:]
    p1_rv_hist  = float.( mapreduce( permutedims, hcat, p1_rv_hist )'  )
    p2_rv_hist  = float.( mapreduce( permutedims, hcat, p2_rv_hist )'  )
    rv_ref_hist = float.( mapreduce( permutedims, hcat, rv_ref_hist )' ) 

    @exfiltrate 

    return tt_hist, p1_rv_hist, p2_rv_hist, rv_ref_hist  
end 

export p_rv_ref_hist 


## ==================================================================== 

function dist_ref_norm( game, params ) 

    _, p1_rv_hist, p2_rv_hist, rv_ref_hist = p_rv_ref_hist( game, params )

    p1_r_diff = p1_rv_hist[:,1:3] - rv_ref_hist[:,1:3] 
    p2_r_diff = p2_rv_hist[:,1:3] - rv_ref_hist[:,1:3] 

    p1_ref_norm = [ norm( p1_r_diff[ii,:] ) for ii in 1 : size(p1_r_diff, 1) ] 
    p2_ref_norm = [ norm( p2_r_diff[ii,:] ) for ii in 1 : size(p2_r_diff, 1) ] 
    
    return p1_ref_norm, p2_ref_norm 
end 

export dist_ref_norm 


## ==================================================================== 

function dist_norm( game, params ) 

    _, p1_rv_hist, p2_rv_hist, _ = p_rv_ref_hist( game, params )

    r_diff = p1_rv_hist[:,1:3] - p2_rv_hist[:,1:3] 
    r_norm = [ norm( r_diff[ii,:] ) for ii in 1 : size(r_diff, 1) ] 
    
    return r_norm
end 

export dist_norm 


## ====================================================================

using JLD2 

function save_games_vec( games_vec, params = games_vec[1].params[1] ) 

    # save_folder = string( "test/results/", params.strategy, "/" ) 
    save_folder = string( "test/results/", "p1_", params.strategy, "_p2_", params.p2_strategy, "/" )  

    if !isdir(save_folder)
        mkdir(save_folder)
    end 

    N_games = length( games_vec ) 

    filename      = string( "games_", N_games, ".jld2" ) 
    full_filename = string( save_folder, filename ) 

    @save full_filename games_vec 

end 

export save_games_vec 


## ====================================================================

function load_games_vec( N_games, p1_strategy = "mixed", p2_strategy = "mixed" ) 

    # save_folder = string( "test/results/", p1_strategy, "/" ) 
    save_folder = string( "test/results/", "p1_", p1_strategy, "_p2_", p2_strategy, "/" )  

    filename      = string( "games_", N_games, ".jld2" ) 
    full_filename = string( save_folder, filename ) 

    @load full_filename games_vec 

    return games_vec 
end 

export load_games_vec 


## ====================================================================

function run_game( rng, N_replan = 10, p1_strategy = "mixed", p2_strategy = "mixed" ) 

    params, players, game = init_game( rng, p1_strategy, p2_strategy )  

    for ii = 1 : N_replan - 1 
        println( "step: ", ii + 1, "\n" ) 
        game = prop_game_step( game, params, rng ) 
    end 

    return game, params 
end 

export run_game 


## ====================================================================

function run_MC_games( rng, N_games, N_replan, p1_strategy = "mixed", p2_strategy = "mixed" ) 

    games_vec = [] 
    
    for jj = 1 : N_games 

        params, players, game = init_game( rng, p1_strategy, p2_strategy ) 

        for ii = 1 : N_replan - 1 
            println( "game: ", jj, " step: ", ii + 1, "\n" ) 
            game = prop_game_step( game, params, rng ) 
        end 

        push!( games_vec, game ) 

    end 
    
    # print_MC_stats( games_vec ) 
    save_games_vec( games_vec )  

    return games_vec 
end 

export run_MC_games 

# ==================================================================== 

using Base.Threads

function run_MC_games_parallel( rng, N_games, N_replan, p1_strategy = "mixed", p2_strategy = "mixed" )

    # Create thread-safe storage
    games_vec = Vector{Any}(undef, N_games)
    
    # Create independent RNGs for each game (thread-safe)
    seeds = rand(rng, UInt, N_games)

    # run games in parallel
    Threads.@threads for i_game = 1:N_games

        local_rng = MersenneTwister(seeds[i_game])
        
        params, players, game = init_game(local_rng, p1_strategy, p2_strategy)

        # propogate each game forward one step at a time 
        for i_step = 1 : N_replan - 1
            println("game: ", i_game, " step: ", i_step + 1)
            game = prop_game_step(game, params, local_rng)
        end

        games_vec[i_game] = game

    end

    save_games_vec(games_vec)

    return games_vec
end

export run_MC_games_parallel 


## ====================================================================

function MC_stats( games_vec, params = games_vec[1].params[1] ) 

    dist_rnorm_all   = [] 
    p1_Unorm_sum_all = [] 
    p2_Unorm_sum_all = [] 
    p1_ref_norm_all  = [] 
    p2_ref_norm_all  = [] 
    
    for ii in eachindex( games_vec ) 

        game = games_vec[ii] 

        # compute norm of U and distance vectors 
        dist_rnorm = dist_norm( game, params ) 
        p1_U_norm, p2_U_norm = U_norm( game, params ) 
        p1_Unorm_sum = cumsum( p1_U_norm ) 
        p2_Unorm_sum = cumsum( p2_U_norm ) 

        p1_ref_norm, p2_ref_norm = dist_ref_norm( game, params ) 

        push!( dist_rnorm_all,      dist_rnorm ) 
        push!( p1_Unorm_sum_all,    p1_Unorm_sum ) 
        push!( p2_Unorm_sum_all,    p2_Unorm_sum ) 
        push!( p1_ref_norm_all,     p1_ref_norm ) 
        push!( p2_ref_norm_all,     p2_ref_norm ) 

    end 

    dist_rnorm_all    = vv2m( dist_rnorm_all ) 
    p1_Unorm_sum_all  = vv2m( p1_Unorm_sum_all ) 
    p2_Unorm_sum_all  = vv2m( p2_Unorm_sum_all ) 
    p1_ref_norm_all   = vv2m( p1_ref_norm_all ) 
    p2_ref_norm_all   = vv2m( p2_ref_norm_all ) 

    dist_norm_mean    = mean( dist_rnorm_all,   dims = 1 )[:] 
    p1_Unorm_sum_mean = mean( p1_Unorm_sum_all, dims = 1 )[:] 
    p2_Unorm_sum_mean = mean( p2_Unorm_sum_all, dims = 1 )[:] 
    p1_ref_norm_mean  = mean( p1_ref_norm_all,  dims = 1 )[:] 
    p2_ref_norm_mean  = mean( p2_ref_norm_all,  dims = 1 )[:] 

    dist_norm_std     = std( dist_rnorm_all,    dims = 1 )[:] 
    p1_Unorm_sum_std  = std( p1_Unorm_sum_all,  dims = 1 )[:] 
    p2_Unorm_sum_std  = std( p2_Unorm_sum_all,  dims = 1 )[:] 
    p1_ref_norm_std   = std( p1_ref_norm_all,   dims = 1 )[:] 
    p2_ref_norm_std   = std( p2_ref_norm_all,   dims = 1 )[:] 

    stats = ( dist_rnorm_all    = dist_rnorm_all, 
              p1_Unorm_sum_all  = p1_Unorm_sum_all, 
              p2_Unorm_sum_all  = p2_Unorm_sum_all, 
              p1_ref_norm_all   = p1_ref_norm_all, 
              p2_ref_norm_all   = p2_ref_norm_all, 
              dist_norm_mean    = dist_norm_mean, 
              p1_Unorm_sum_mean = p1_Unorm_sum_mean, 
              p2_Unorm_sum_mean = p2_Unorm_sum_mean, 
              p1_ref_norm_mean  = p1_ref_norm_mean, 
              p2_ref_norm_mean  = p2_ref_norm_mean, 
              dist_norm_std     = dist_norm_std, 
              p1_Unorm_sum_std  = p1_Unorm_sum_std, 
              p2_Unorm_sum_std  = p2_Unorm_sum_std, 
              p1_ref_norm_std   = p1_ref_norm_std, 
              p2_ref_norm_std   = p2_ref_norm_std )  

    return stats  
end 

export MC_stats 


## ====================================================================

function print_MC_stats( 
    games_vec, 
    params      = games_vec[1].params[1], 
    print_stats = true 
) 

    stats = MC_stats( games_vec, params ) 
    
    dist_norm_mean_mean   = @sprintf "%.3g" mean(stats.dist_norm_mean)  
    p1_Unorm_mean_end     = @sprintf "%.3g" stats.p1_Unorm_sum_mean[end]  
    p2_Unorm_mean_end     = @sprintf "%.3g" stats.p2_Unorm_sum_mean[end]  
    p1_ref_norm_mean_mean = @sprintf "%.3g" mean(stats.p1_ref_norm_mean)
    p2_ref_norm_mean_mean = @sprintf "%.3g" mean(stats.p2_ref_norm_mean) 

    sprintf_stats = ( dist_norm_mean_mean = dist_norm_mean_mean, 
                      p1_Unorm_mean_end   = p1_Unorm_mean_end, 
                      p2_Unorm_mean_end   = p2_Unorm_mean_end, 
                      p1_ref_norm_mean_mean = p1_ref_norm_mean_mean, 
                      p2_ref_norm_mean_mean = p2_ref_norm_mean_mean )  

    if print_stats == true 

        println( "games = ", length(games_vec) )
        println( "p1 strategy = ", params.strategy, ", p2 strategy = ", params.p2_strategy ) 
        println( "mean player distance: ", dist_norm_mean_mean ) 
        println( "mean cumsum norm of U vectors: p1 = ", p1_Unorm_mean_end, ", p2 = ", p2_Unorm_mean_end ) 
        println( "mean player distance from reference orbit: p1 = ", p1_ref_norm_mean_mean, ", p2 = ", p2_ref_norm_mean_mean )      

    end 

    return sprintf_stats 
end 

export print_MC_stats 







