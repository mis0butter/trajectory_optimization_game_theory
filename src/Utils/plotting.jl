## ============================================ ##
# plot Cartesian axes 

"Plot x, y, and z Cartesian axes using GLMakie "
function plot_axes3d( 
    r   = 6378.0 / 3,   # radius of axes 
    fig = nothing,      # figure handle 
) 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    xyz = [ zeros(3) for i in 1:3 ] 
    uvw = r .* [ [1,0,0] , [0,1,0] , [0,0,1] ] 

    width = r/50 
    fig = plot_vector3d( [ xyz[1] ] , [ uvw[1] ], nothing, width, :red ) 
    fig = plot_vector3d( [ xyz[2] ] , [ uvw[2] ], fig, width, :blue ) 
    fig = plot_vector3d( [ xyz[3] ] , [ uvw[3] ], fig, width, :green  )  

    return fig 
end

export plot_axes3d 

## ============================================ ## 

""" 
Plot a 3D line using GLMakie 

Example usage: 

    x = collect( range(-pi, pi, 100) ) 
    y = sin.(x) 
    z = cos.(x) 

    fig = plot_3d( x, y, z )
"""

function plot_line3d( xyz, fig = nothing ) 
    plot_line3d( xyz[:,1], xyz[:,2], xyz[:,3], fig ) 
end 

function plot_line3d( 
    x,              # [N,1] grid of points 
    y,              # [N,1] grid of points 
    z,              # [N,1] grid of points  
    fig = nothing,  # figure handle 
    color = :black   # line color 
) 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    # plot orbit 
    lines!( x, y, z; linewidth = 2, alpha = 0.5, color = color  ) 

    return fig 
end 
    
export plot_line3d  

## ============================================ ## 

""" 
Plot an orbit using GLMakie 

Example usage: 

    x = collect( range(-pi, pi, 100) ) 
    y = sin.(x) 
    z = cos.(x) 

    fig = plot_orbit( [x y z] )
"""

function plot_orbit( 
    rv,                 # [N,3] matrix of state vectors 
    fig    = nothing,   # figure handle 
    labels = false      # boolean for labeling start and end points 
) 

    text_offset = (0,10) 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    # if isnothing(fig) 
    #     fig = Figure() 
    #     Axis3(fig[1, 1], aspect=DataAspect(), 
    #         xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)", 
    #         title = "Transfer Solution") 
    # end 

    # plot orbit 
    lines!( rv[:,1], rv[:,2], rv[:,3]; linewidth = 2 ) 
    scatter!( rv[1,1], rv[1,2], rv[1,3]; marker = :circle, markersize = 10, color = :black ) 
    scatter!( rv[end,1], rv[end,2], rv[end,3]; marker = :utriangle, markersize = 10, color = :black ) 

    # add labels 
    if labels 
        text!( rv[1,1], rv[1,2], rv[1,3]; text = "start", color = :gray, offset = text_offset, align = (:center, :bottom) ) 
        text!( rv[end,1], rv[end,2], rv[end,3]; text = "end", color = :gray, offset = text_offset, align = (:center, :bottom) ) 
    end 

    Auto() 

    return fig 
end 
    
export plot_orbit 

## ============================================ ## 

# colormap options: 
#   jblue 
#   copper 
#   diverging_tritanopic_cwr_75_98_c20_n256 <-- this one 

"""
Plot a surface with a colorbar using GLMakie 

Example usage: 

    x = y = range(-pi, pi, 100)
    z = sin.(x) .* cos.(y') 

    fig = plot_surface( x, y, z ) 
"""

function plot_surface( 
    x,                  # [N,1] grid of points 
    y,                  # [N,1] grid of points 
    z,                  # [N,N] grid of points evaluated at x and y 
    fig   = nothing,    # figure handle 
    alpha = 1.0,        # transparency 
) 

    fignothing = false 
    if isnothing(fig) 
        fignothing = true 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    cmap = ( :diverging_tritanopic_cwr_75_98_c20_n256, alpha )
    hm   = GLMakie.surface!( x, y, z, colormap = cmap ) 

    if fignothing 
        Colorbar( fig[1,2], hm, height = Relative(0.5) )
    end 

    return fig 
end 

export plot_surface 

## ============================================ ##

"""
Plot scatter using GLMakie

Example usage: 

    x = y = range(-pi, pi, 100)
    z = sin.(x) .* cos.(y') 

    fig = plot_surface( x, y, z ) 
    fig = plot_contour3d( x, y, z ) 
    fig = plot_scatter3d( x, y, z ) 
"""

function plot_scatter3d( xyz, fig = nothing ) 
    plot_scatter3d( xyz[:,1], xyz[:,2], xyz[:,3], fig ) 
end 

function plot_scatter3d( 
    x,                      # [N,1] grid of points 
    y,                      # [N,1] grid of points 
    z,                      # [N,N] grid of points evaluated at x and y 
    fig    = nothing,       # figure handle 
    marker = :utriangle,    # marker type 
    color  = :black,        # marker color 
    markersize = 12,        # marker size 
    text   = nothing,       # text to add to plot 
) 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    if isequal(length(z), 1)
        GLMakie.scatter!( x, y, z, marker = marker, markersize = markersize, color = color, strokecolor = color ) 
        if !isnothing(text) 
            text!( x, y, z; text = text, color = :black, offset = (0,15), align = (:center, :bottom) ) 
        end
    else 
        hm = GLMakie.scatter!( x, y, z, markersize = 5, color = color, strokecolor = color ) 
    end 

    return fig 
end 

export plot_scatter3d 

## ============================================ ##

""" 
Plot a contour with a colorbar using GLMakie

Example usage: 

    x = y = range(-pi, pi, 100)
    z = sin.(x) .* cos.(y') 

    fig = plot_contour3d( x, y, z ) 
""" 

function plot_contour3d( 
    x,              # [N,1] grid of points 
    y,              # [N,1] grid of points 
    z,              # [N,N] grid of points evaluated at x and y 
    fig = nothing,  # figure handle 
    levels = 20,    # number of contour levels 
) 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    hm  = GLMakie.contour3d!(x, y, z, levels = levels) 

    if fignothing 
        clim = ( minimum(z), maximum(z) ) 
        Colorbar( fig[1, 2], limits = clim, height = Relative(0.5) )
    end 

    return fig 
end 

export plot_contour3d 

## ============================================ ##

"""
Plot vector using GLMakie. 

Example usage: 

    r   = 6378.0
    xyz = [ zeros(3) for i in 1:3 ] 
    uvw = r .* [ [1,0,0] , [0,1,0] , [0,0,1] ] 

    fig = plot_vector3d( [ xyz[1] ] , [ uvw[1] ], nothing, r/100, :red ) 
    fig = plot_vector3d( [ xyz[2] ] , [ uvw[2] ], fig, r/100, :blue ) 
    fig = plot_vector3d( [ xyz[3] ] , [ uvw[3] ], fig, r/100, :green ) 
"""

"""
    plot_vector3d(origin, direction; fig, width, color, text)

Plot 3D vector arrow(s) from origin point(s) in given direction(s).

# Arguments
- `origin`: Single point [x,y,z] or multiple points [[x1,y1,z1], [x2,y2,z2], ...]
- `direction`: Single vector [u,v,w] or multiple vectors [[u1,v1,w1], ...]
- `fig`: Existing figure handle (creates new if nothing)
- `width`: Arrow line width 
- `color`: Arrow color
- `text`: Optional label (only for single vector)
"""
function plot_vector3d( 
    origin,                     
    direction,                  
    fig    = nothing,           
    width  = nothing,           
    color  = :black,            
    text   = nothing,           
) 
    # normalize inputs to vectors of SVectors
    origins    = _to_point_list(origin)
    directions = _to_point_list(direction)
    
    # default width based on direction magnitude
    if isnothing(width)
        width = norm(directions[1]) / 50
    end

    # convert to Makie types
    ps = [Point3f(p...) for p in origins]
    ns = [Vec3f(d...)   for d in directions]

    # create figure if needed
    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1], aspect = :data) 
    end 

    arrows!(  
        ps, ns, 
        fxaa = true,
        linecolor = color, arrowcolor = color,
        linewidth = width, arrowsize = Vec3f(3*width, 3*width, 4*width),
        align = :origin, 
    )

    if !isnothing(text) 
        length(origins) == 1 || error("Can only label one vector at a time.")
        tip = origins[1] .+ directions[1]
        text!(tip[1], tip[2], tip[3]; text = text, color = :black, offset = (0, 15), align = (:center, :bottom)) 
    end

    return fig 
end

# Helper: convert various input formats to Vector of SVectors
function _to_point_list(input)
    # single point as tuple or SVector: (x,y,z) or SVector(x,y,z)
    if input isa Tuple || input isa StaticVector
        return [SVector{3}(input...)]
    end
    
    # single point as 1D vector: [x, y, z]
    if input isa AbstractVector && eltype(input) <: Number
        return [SVector{3}(input...)]
    end
    
    # matrix: each row is a point (N×3)
    if input isa AbstractMatrix
        return [SVector{3}(input[i, :]...) for i in 1:size(input, 1)]
    end
    
    # already a list of points: [[x1,y1,z1], [x2,y2,z2], ...]
    if input isa AbstractVector
        return [SVector{3}(p...) for p in input]
    end
    
    error("Unsupported input format for plot_vector3d")
end 

export plot_vector3d 

## ============================================ ##

"Plot propagated orbit with delta v using GLMakie "
function plot_prop_Δv(  
    rv_0,               # initial state vector 
    Δv_sol,             # [N,3] Δv vector 
    N,                  # number of segments 
    tof_N_sol,          # time of flight for each segment 
    mu  = 1.0,          # gravitational parameter 
    fig = nothing,      # figure handle 
    plot_vector = true, # plot Δv vector 
)

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    # propagate 2 body 
    t, rv_2Body = prop_2Body_tof_Nseg( rv_0, Δv_sol, N, tof_N_sol, mu ) 

    # propagate kepler 
    t, rv_kepler = prop_kepler_tof_Nseg( rv_0, Δv_sol, N, tof_N_sol, mu ) 

    # plot 
    fig = plot_orbit( rv_2Body, fig ) 
    # fig = plot_vector3d( [ x0_P[1:3] ], 500 * [ Δv ], fig ) 

    if plot_vector 

        # set up vector plotting 
        nodes_N = rv_kepler[1:N, 1:3] 
        xyz     = copy(nodes_N) 
        uvw     = copy(2000 * Δv_sol)

        fig = plot_vector3d( xyz, uvw, fig, 100 ) 

    end

    return fig 
end 

export plot_prop_Δv 

## ============================================ ##

"Plot lines of polygon at rv input"
function plot_polygon( 
    rv_state,                   # [N,6] state vector 
    parameters,                 # struct of parameters 
    figure = plot_axes3d()      # figure handle 
    ) 

    R_polygon = parameters.R_polygon 

    vertices  = polygon_vertices( rv_state, parameters ) 
    
    figure = plot_scatter3d( rv_state[1], rv_state[2], rv_state[3], figure ) 
    
    # center of polygon 
    r_vec = rv_state[1:3] 
    axis_1, axis_2, axis_3 = axis_123( rv_state ) 

    # ok, let's plot this so that it all looks right 
    # fig = plot_vector3d( [ r_vec ] , [ axis_1 * r ] , fig, r/100, :black, "1" ) 
    # fig = plot_vector3d( [ r_vec ] , [ axis_2 * r ] , fig, r/100, :black, "2" ) 
    # fig = plot_vector3d( [ r_vec ] , [ axis_3 * r ] , fig, r/100, :black, "3" ) 

    # vertices of polygon along axis 2-3 plane 
    # ok, let's define the distance of vertices of polygon from center: how about r / 100 ? 

    # top vertex: move up from r_f along axis 3 
    r_top = r_vec + axis_3 * R_polygon 
    # fig   = plot_scatter3d( r_top[1], r_top[2], r_top[3], fig, :circle ) 

    # top-inner vertex: move up from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = cosd(60) * axis_3 * R_polygon + sind(60) * axis_2 * R_polygon
    r_topin  = r_vec + vec
    # fig      = plot_scatter3d( r_topin[1], r_topin[2], r_topin[3], fig, :circle )  

    mat = [ r_top' ; r_topin' ] 
    figure = plot_line3d( mat, figure ) 

    # bottom-inner vertex: move down from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * R_polygon + sind(60) * axis_2 * R_polygon
    r_botin  = r_vec + vec 
    # fig      = plot_scatter3d( r_botin[1], r_botin[2], r_botin[3], fig, :circle ) 

    mat = [ r_topin' ; r_botin' ] 
    figure = plot_line3d( mat, figure ) 

    # bottom vertex: move down from r_f along axis 3 
    r_bot = r_vec - axis_3 * R_polygon
    # fig   = plot_scatter3d( r_bot[1], r_bot[2], r_bot[3], fig, :circle ) 

    mat = [ r_botin' ; r_bot' ] 
    figure = plot_line3d( mat, figure ) 

    # bottom-outer vertex: move down from r_f along axis 3 and right along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * R_polygon - sind(60) * axis_2 * R_polygon
    r_botout = r_vec + vec 
    # fig      = plot_scatter3d( r_botout[1], r_botout[2], r_botout[3], fig, :circle ) 

    mat = [ r_bot' ; r_botout' ] 
    figure = plot_line3d( mat, figure ) 

    # top-outer vertex: move up from r_f along axis 3 and right along axis 2, 60 degrees 
    vec       = cosd(60) * axis_3 * R_polygon - sind(60) * axis_2 * R_polygon
    r_topout  = r_vec + vec 
    # fig       = plot_scatter3d( r_topout[1], r_topout[2], r_topout[3], fig, :circle ) 

    mat = [ r_botout' ; r_topout' ] 
    figure = plot_line3d( mat, figure ) 

    mat = [ r_topout' ; r_top' ] 
    figure = plot_line3d( mat, figure ) 

    return figure 
end 

export plot_polygon 

## ============================================ ## 

""" 
Plot a candidate trajectory using GLMakie 

Example usage: 

    x = collect( range(-pi, pi, 100) ) 
    y = sin.(x) 
    z = cos.(x) 

    fig = plot_traj_cand( [x y z] )
"""

function plot_traj_cand( 
    rv,                 # [N,3] matrix of state vectors 
    color  = 1,         # color  
    alpha  = 1,         # transparency 
    linestyle = :dash,  # line style 
    linewidth = 2,      # line width 
    fig    = nothing,   # figure handle 
) 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    # plot orbit 
    lines!( rv[:,1], rv[:,2], rv[:,3]; linewidth = linewidth, color = color, alpha = alpha, linestyle = linestyle ) 
    scatter!( rv[1,1], rv[1,2], rv[1,3]; marker = :circle, markersize = 20, color = color ) 
    # scatter!( rv[end,1], rv[end,2], rv[end,3]; marker = :utriangle, markersize = 10, color = :black ) 

    Auto() 

    return fig 
end 
    
export plot_traj_cand 

## ============================================ ##

"Plot propagated orbit with delta v using GLMakie "
function plot_Δv_weights(  
    game,               # game struct 
    params,             # struct of parameters 
    k,                  # k_replan step of the game      
    fig = nothing,      # figure handle 
)

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    n_vertices = size( game.p1_state[1].cost , 1 ) 

    # get params 
    N  = params.N 
    mu = params.mu 
    tof_N_sol = params.tof / params.N 

    p1_chosen = game.p1_state[k].chosen ; p1_color = :blue 
    p2_chosen = game.p2_state[k].chosen ; p2_color = :red 

    for i in 1 : n_vertices  
        
        # ----------------------- #
        # propagate 2 body and plot 

        rv_0   = game.p1_state[k].rv_0_hist[1,:] 
        Δv_sol = game.p1_state[k].U[i] 
        weight = game.p1_state[k].weights[i] 

        t, rv_2Body = prop_2Body_tof_Nseg( rv_0, Δv_sol, N, tof_N_sol, mu ) 
        fig = plot_traj_cand( rv_2Body, p1_color, weight, :dash, 5, fig ) 
        if i == p1_chosen 
            fig = plot_traj_cand( rv_2Body, p1_color, weight, :dot, 2, fig ) 
        end 
        
        # ----------------------- #
        # propagate 2 body and plot 

        rv_0   = game.p2_state[k].rv_0_hist[1,:] 
        Δv_sol = game.p2_state[k].U[i] 
        weight = game.p2_state[k].weights[i] 

        t, rv_2Body = prop_2Body_tof_Nseg( rv_0, Δv_sol, N, tof_N_sol, mu ) 
        fig = plot_traj_cand( rv_2Body, p2_color, weight, :dash, 5, fig ) 
        if i == p2_chosen 
            fig = plot_traj_cand( rv_2Body, p2_color, weight, :dot, 2, fig ) 
        end 

    end 

    return fig 
end 

export plot_Δv_weights 

## ============================================ ##

function plot_p1_p2_traj( game, parameters, kk ) 

    # propagate SC state forward 
    p1 = game.p1_state[ kk ] 
    p2 = game.p2_state[ kk ] 

    # why tf was I doing this 
    # player1_chosen, player2_chosen = p_strategy( game, kk, parameters.strategy ) 

    p1_chosen = p1.chosen 
    p2_chosen = p2.chosen 

    # get current state 
    rv_E = p1.X[ p1_chosen ][ parameters.k_tt_replan + 1, : ]
    rv_P = p2.X[ p2_chosen ][ parameters.k_tt_replan + 1, : ]

    # plot 
    fig = plot_axes3d(  ) 
    # fig = plot_orbit( rv_E_hist, fig ) 
    # fig = plot_orbit( rv_P_hist, fig ) 
    fig = plot_polygon( game.rv_ref_E[kk][end,:], parameters, fig ) 

    for i in 1 : parameters.k_tt_replan 

        # plot player 1 
        rv_E = p1.X[ p1.chosen ][ i, : ] 
        fig  = plot_scatter3d( rv_E[1], rv_E[2], rv_E[3], fig, :circle, :blue, 20 ) 

        # plot player 2 
        rv_P = p2.X[ p2.chosen ][ i, : ] 
        fig  = plot_scatter3d( rv_P[1], rv_P[2], rv_P[3], fig, :circle, :red, 20 ) 

    end 

    # plot triangles on chosen vertices 
    rv_E_target = p1.X[ p1.chosen ][ end, : ] 
    rv_P_target = p2.X[ p2.chosen ][ end, : ] 
    fig = plot_scatter3d( rv_E_target[1], rv_E_target[2], rv_E_target[3], fig, :utriangle, :blue, 15 ) 
    fig = plot_scatter3d( rv_P_target[1], rv_P_target[2], rv_P_target[3], fig, :utriangle, :red, 15 ) 
    
    # if k > 1 
    if kk > 1 
        for j = 1 : kk - 1 

            p1 = game.p1_state[ j ] 
            p2 = game.p2_state[ j ] 

            rv_E = p1.X[ p1.chosen ][ 1 : parameters.k_tt_replan + 1, : ] 
            lines!( rv_E[:,1], rv_E[:,2], rv_E[:,3]; linewidth = 2, color = :blue ) 

            rv_P = p2.X[ p2.chosen ][ 1 : parameters.k_tt_replan + 1, : ] 
            lines!( rv_P[:,1], rv_P[:,2], rv_P[:,3]; linewidth = 2, color = :red ) 

        end 
    end 
    
    # plot all the weights 
    fig = plot_Δv_weights( game, parameters, kk, fig ) 

    # title 
    ax = fig.current_axis 

    # title_string = string( parameters.strategy, " game: k_replan = ", kk ) 
    title_string = string( "p1 strategy = ", parameters.strategy, ", p2 strategy = ", parameters.p2_strategy, "\n k_replan = ", kk ) 
    ax.x.title = title_string 

    return fig 
end 

export plot_p1_p2_traj 


## ============================================ ##

function plot_games_costs( fig, x_fig, y_fig, games_vec ) 

    # get costs 
    games_costs, games_costs_mean, games_costs_std = stage_cost_games_fn( games_vec )

    # mean and mean-std strings 
    string_mean = @sprintf "%.3g" mean(games_costs_mean) 
    string_std_mean  = @sprintf "%.3g" mean(games_costs_std) 

    # get time vector 
    tt_hist, _, _, _ = p_rv_ref_hist( games_vec[1] ) 

    # mean +/- std 
    y_upper = games_costs_mean .+ games_costs_std 
    y_lower = games_costs_mean .- games_costs_std 

    # axis title 
    title_string = "stage costs (game value)"
    title_string = string( "mean cost = ", string_mean, ", mean std = ", string_std_mean )  

    # create axis 
    ax = Axis( fig[x_fig, y_fig], xlabel = "time (s)", title = title_string )

    # plot mean 
    lines!( ax, tt_hist[1:end-1], games_costs_mean, color = :green )   

    # plot mean +/- std ribbons 
    fill_between!(ax, tt_hist[1:end-1], y_lower, y_upper, color = :green, alpha = 0.25 ) 

    # plot individual games 
    for ii in eachindex(games_vec)
        lines!( ax, tt_hist[1:end-1], games_costs[ii,:][:], color = :green, alpha = 0.1 ) 
    end 

    return fig 
end 


## ============================================ ##

function plot_ref_stats( fig, x_fig, y_fig, games_vec ) 

    # get stats 
    stats = MC_stats( games_vec ) 
    sprintf_stats = print_MC_stats( games_vec ) 

    # get time vector 
    tt, _, _, _ = p_rv_ref_hist( games_vec[1] ) 

    # temp strings for title string  
    temp1 = @sprintf "%.3g" mean(stats.p1_ref_norm_std) 
    temp2 = @sprintf "%.3g" mean(stats.p2_ref_norm_std) 

    # axis title string 
    title_string = string( 
        "mean player distance from ref orbit: ", 
        "\n p1 mean = ",  sprintf_stats.p1_ref_norm_mean_mean, 
        ", std = ", temp1, 
        "\n p2 mean = ",  sprintf_stats.p2_ref_norm_mean_mean, 
        ", std = ", temp2  
    ) 

    # create axis 
    ax3 = Axis( fig[x_fig, y_fig], xlabel = "time", title = title_string ) 

    # plot reference norm mean 
    p1_ax3 = lines!( ax3, tt, stats.p1_ref_norm_mean, color = :blue ) 
    p2_ax3 = lines!( ax3, tt, stats.p2_ref_norm_mean, color = :red ) 

    # plot individual games 
    for ii in eachindex(games_vec)
        lines!( ax3, tt, stats.p1_ref_norm_all[ii,:][:], color = :blue, alpha = 0.1 ) 
        lines!( ax3, tt, stats.p2_ref_norm_all[ii,:][:], color = :red, alpha = 0.1 ) 
    end 

    return fig 
end 


## ============================================ ##

function plot_player_distance( fig, x_fig, y_fig, games_vec ) 

    # get stats 
    stats = MC_stats( games_vec ) 
    sprintf_stats = print_MC_stats( games_vec ) 

    # get time vector and params 
    params = games_vec[1].params[1] 
    tt, _, _, _ = p_rv_ref_hist( games_vec[1] ) 
    
    # title string 
    title_string = string( 
        "p1 strategy = ",   params.strategy, 
        ", p2 strategy = ", params.p2_strategy, 
        "\n mean player distance = ", sprintf_stats.dist_norm_mean_mean, 
        ", mean std = ", @sprintf "%.3g" mean(stats.dist_norm_std) 
    )

    # create axis 
    ax2 = Axis( fig[x_fig, y_fig], xlabel = "time (s)", title = title_string )  

    # plot mean +/- std ribbons 
    y_upper = stats.dist_norm_mean .+ stats.dist_norm_std 
    y_lower = stats.dist_norm_mean .- stats.dist_norm_std 
    fill_between!(ax2, tt, y_lower, y_upper, color = :green, alpha = 0.25 )

    # plot mean 
    lines!( ax2, tt, stats.dist_norm_mean, color = :green ) 

    # plot individual games 
    for ii in eachindex(games_vec)
        lines!( ax2, tt, stats.dist_rnorm_all[ii,:][:], color = :green, alpha = 0.1 ) 
    end 

    return fig 
end 


## ============================================ ##

function plot_cumsum_U( fig, x_fig, y_fig, games_vec ) 

    # get stats 
    stats = MC_stats( games_vec ) 
    sprintf_stats = print_MC_stats( games_vec ) 

    # get time vector and params 
    params = games_vec[1].params[1] 
    tt, _, _, _ = p_rv_ref_hist( games_vec[1] ) 
    
    # title string 
    title_string = string( 
        "mean cumsum norm of U vectors \n", 
        "p1 = ",    sprintf_stats.p1_Unorm_mean_end, 
        ", p2 = ",  sprintf_stats.p2_Unorm_mean_end 
    )  

    # create axis 
    ax1 = Axis( fig[x_fig, y_fig], xlabel = "time (s)", title = title_string ) 

    # plot cumsum norm of U vectors 
    p1_ax1 = lines!( ax1, tt[ 1 : end - 1 ], stats.p1_Unorm_sum_mean, color = :blue ) 
    p2_ax1 = lines!( ax1, tt[ 1 : end - 1 ], stats.p2_Unorm_sum_mean, color = :red ) 

    # legend 
    Legend( fig[x_fig, y_fig + 1], [ p1_ax1, p2_ax1 ], ["p1", "p2"] ) 

    # plot individual games 
    for ii in eachindex(games_vec)
        lines!( ax1, tt[1:end-1], stats.p1_Unorm_sum_all[ii,:][:], color = :blue, alpha = 0.1 ) 
        lines!( ax1, tt[1:end-1], stats.p2_Unorm_sum_all[ii,:][:], color = :red,  alpha = 0.1 ) 
    end 

    return fig 
end 


## ============================================ ##

function plot_MC_stats( games_vec, params = games_vec[1].params[1] )

    # figure 
    fig = Figure( size = (1000, 600) ) 

    # top left 
    fig = plot_player_distance( fig, 1, 1, games_vec )  

    # top right 
    fig = plot_cumsum_U( fig, 1, 2, games_vec ) 

    # bottom left 
    fig = plot_games_costs( fig, 2, 1, games_vec ) 

    # bottom right 
    fig = plot_ref_stats( fig, 2, 2, games_vec ) 

    return fig 
end 

export plot_MC_stats 

## ============================================ ## 

function plot_game_stats( game, params ) 

    r_norm = dist_norm( game, params ) 
    p1_U_norm, p2_U_norm = U_norm( game, params ) 
    p1_Unorm_sum = cumsum( p1_U_norm ) 
    p2_Unorm_sum = cumsum( p2_U_norm ) 

    fig = Figure( size = (600, 600) )

    # title_string = string( params.strategy, " game \n norm of U vectors" ) 
    title_string = string( "p1 strategy = ", params.strategy, ", p2 strategy = ", params.p2_strategy, "\n norm of U vectors" ) 

    ax1 = Axis( fig[1,1], xlabel = "time", title = title_string ) 
    p1_ax1 = lines!( ax1, 1 : length(p1_U_norm), p1_U_norm, color = :blue ) 
    p2_ax1 = lines!( ax1, 1 : length(p2_U_norm), p2_U_norm, color = :red ) 
    Legend( fig[1,2], [ p1_ax1, p2_ax1 ], ["p1", "p2"] ) 

    p1_Unorm_end = @sprintf "%.3g" p1_Unorm_sum[end] 
    p2_Unorm_end = @sprintf "%.3g" p2_Unorm_sum[end] 
    title_string = string( "cumsum of U norm \n p1 = ", p1_Unorm_end, ", p2 = ", p2_Unorm_end )  
    ax2 = Axis( fig[2,1], xlabel = "time", title = title_string ) 
    lines!( ax2, 1 : length(p1_Unorm_sum), p1_Unorm_sum, color = :blue ) 
    lines!( ax2, 1 : length(p2_Unorm_sum), p2_Unorm_sum, color = :red ) 

    ax3 = Axis( fig[3,1], xlabel = "time", title = "player distance" ) 
    lines!( ax3, 1 : length(r_norm), r_norm, color = :green )  

    @exfiltrate 

    return fig 
end 

export plot_game_stats 


