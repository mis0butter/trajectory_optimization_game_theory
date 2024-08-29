using GLMakie 

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

function plot_vector3d( 
    xyz,                        # [N] vector of (x,y,z) origin points 
    uvw,                        # [N] vector of (u,v,w) vector directions 
    fig    = nothing,           # figure handle 
    width  = norm(uvw[1])/100,  # arrow width 
    color  = :black,            # marker color 
    text   = nothing,           # text to add to plot 
) 

    # check type --> must be vectors of vectors 
    if xyz isa AbstractMatrix 
        xyz = m2vv(xyz)
    end 
    if uvw isa AbstractMatrix 
        uvw = m2vv(uvw)
    end 

    # adjust because stupid arrows plots the tails at the Point3f points 
    xyz += uvw 

    # convert to Points3f and Vec3f for arrows function 
    ps  = [ Point3f(x,y,z) for (x,y,z) in xyz ] 
    ns  = [ Vec3f(x,y,z) for (x,y,z) in uvw ] 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1]) 
    end 

    arrows!(  
        ps, ns, fxaa = true, # turn on anti-aliasing
        linecolor = color, arrowcolor = color,
        linewidth = width, arrowsize = 2 * width .* Vec3f(1, 1, 1),
        align = :center, 
    )

    if !isnothing(text) 
        if size(xyz, 1) > 1 
            error("Can only label one vector at time.")
        end 
        x = xyz[1][1] ; y = xyz[1][2] ; z = xyz[1][3] 
        text!( x, y, z; text = text, color = :black, offset = (0,15), align = (:center, :bottom) ) 
    end

    return fig 
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
    rv_vec,                 # [N,6] state vector 
    fig = plot_axes3d(),    # figure handle 
    r   = R_polygon,            # radius of Earth (km) 
) 

    vertices = polygon_vertices( rv_vec ) 
    
    fig = plot_scatter3d( rv_vec[1], rv_vec[2], rv_vec[3], fig ) 
    
    # center of polygon 
    r_vec = rv_vec[1:3] 
    axis_1, axis_2, axis_3 = axis_123( rv_vec ) 

    # ok, let's plot this so that it all looks right 
    # fig = plot_vector3d( [ r_vec ] , [ axis_1 * r ] , fig, r/100, :black, "1" ) 
    # fig = plot_vector3d( [ r_vec ] , [ axis_2 * r ] , fig, r/100, :black, "2" ) 
    # fig = plot_vector3d( [ r_vec ] , [ axis_3 * r ] , fig, r/100, :black, "3" ) 

    # vertices of polygon along axis 2-3 plane 
    # ok, let's define the distance of vertices of polygon from center: how about r / 100 ? 

    # top vertex: move up from r_f along axis 3 
    r_top = r_vec + axis_3 * r 
    # fig   = plot_scatter3d( r_top[1], r_top[2], r_top[3], fig, :circle ) 

    # top-inner vertex: move up from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = cosd(60) * axis_3 * r + sind(60) * axis_2 * r
    r_topin  = r_vec + vec
    # fig      = plot_scatter3d( r_topin[1], r_topin[2], r_topin[3], fig, :circle )  

    mat = [ r_top' ; r_topin' ] 
    fig = plot_line3d( mat, fig ) 

    # bottom-inner vertex: move down from r_f along axis 3 and left along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * r + sind(60) * axis_2 * r
    r_botin  = r_vec + vec 
    # fig      = plot_scatter3d( r_botin[1], r_botin[2], r_botin[3], fig, :circle ) 

    mat = [ r_topin' ; r_botin' ] 
    fig = plot_line3d( mat, fig ) 

    # bottom vertex: move down from r_f along axis 3 
    r_bot = r_vec - axis_3 * r
    # fig   = plot_scatter3d( r_bot[1], r_bot[2], r_bot[3], fig, :circle ) 

    mat = [ r_botin' ; r_bot' ] 
    fig = plot_line3d( mat, fig ) 

    # bottom-outer vertex: move down from r_f along axis 3 and right along axis 2, 60 degrees 
    vec      = - cosd(60) * axis_3 * r - sind(60) * axis_2 * r
    r_botout = r_vec + vec 
    # fig      = plot_scatter3d( r_botout[1], r_botout[2], r_botout[3], fig, :circle ) 

    mat = [ r_bot' ; r_botout' ] 
    fig = plot_line3d( mat, fig ) 

    # top-outer vertex: move up from r_f along axis 3 and right along axis 2, 60 degrees 
    vec       = cosd(60) * axis_3 * r - sind(60) * axis_2 * r
    r_topout  = r_vec + vec 
    # fig       = plot_scatter3d( r_topout[1], r_topout[2], r_topout[3], fig, :circle ) 

    mat = [ r_botout' ; r_topout' ] 
    fig = plot_line3d( mat, fig ) 

    mat = [ r_topout' ; r_top' ] 
    fig = plot_line3d( mat, fig ) 

    return fig 
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

function plot_p1_p2_traj( game, params, kk, strategy = "mixed" ) 

    # propagate SC state forward 
    p1 = game.p1_state[ kk ] 
    p2 = game.p2_state[ kk ] 

    p1_chosen, p2_chosen = p_strategy( game, kk, strategy ) 

    # get current state 
    rv_E = p1.X[ p1_chosen ][ params.k_tt_replan + 1, : ]
    rv_P = p2.X[ p2_chosen ][ params.k_tt_replan + 1, : ]

    # plot 
    fig = plot_axes3d(  ) 
    # fig = plot_orbit( rv_E_hist, fig ) 
    # fig = plot_orbit( rv_P_hist, fig ) 
    fig = plot_polygon( game.rv_ref_E_polygon[kk][end,:], fig ) 

    for i in 1 : params.k_tt_replan 

        # plot player 1 
        rv_E = p1.X[ p1.chosen ][ i, : ] 
        fig = plot_scatter3d( rv_E[1], rv_E[2], rv_E[3], fig, :circle, :blue, 20 ) 

        # plot player 2 
        rv_P = p2.X[ p2.chosen ][ i, : ] 
        fig = plot_scatter3d( rv_P[1], rv_P[2], rv_P[3], fig, :circle, :red, 20 ) 

    end 

    # plot triangles on chosen vertices 
    rv_E_des = p1.X[ p1.chosen ][ end, : ] 
    rv_P_des = p2.X[ p2.chosen ][ end, : ] 
    fig = plot_scatter3d( rv_E_des[1], rv_E_des[2], rv_E_des[3], fig, :utriangle, :blue, 15 ) 
    fig = plot_scatter3d( rv_P_des[1], rv_P_des[2], rv_P_des[3], fig, :utriangle, :red, 15 ) 
    
    # if k > 1 
    if kk > 1 
        for j = 1 : kk - 1 

            p1 = game.p1_state[ j ] 
            p2 = game.p2_state[ j ] 

            rv_E = p1.X[ p1.chosen ][ 1 : params.k_tt_replan + 1, : ] 
            lines!( rv_E[:,1], rv_E[:,2], rv_E[:,3]; linewidth = 2, color = :blue ) 


            rv_P = p2.X[ p2.chosen ][ 1 : params.k_tt_replan + 1, : ] 
            lines!( rv_P[:,1], rv_P[:,2], rv_P[:,3]; linewidth = 2, color = :red ) 

        end 
    end 
    
    # plot all the weights 
    fig = plot_Δv_weights( game, params, kk, fig )      

    # title 
    ax = fig.current_axis 
    ax.x.title = string( "k_replan = ", kk )

    return fig 
end 

export plot_p1_p2_traj 

## ============================================ ##

function plot_MC_stats( gameS, params )

    dist_rnorm_all   = [] 
    p1_Unorm_sum_all = [] 
    p2_Unorm_sum_all = [] 
    for ii in eachindex(gameS)

        game = gameS[ii] 

        # compute norm of U and distance vectors 
        dist_rnorm = dist_norm( game, params ) 
        p1_U_norm, p2_U_norm = U_norm( game, params ) 
        p1_Unorm_sum = cumsum( p1_U_norm ) 
        p2_Unorm_sum = cumsum( p2_U_norm ) 

        push!( dist_rnorm_all, dist_rnorm ) 
        push!( p1_Unorm_sum_all, p1_Unorm_sum ) 
        push!( p2_Unorm_sum_all, p2_Unorm_sum ) 

    end 

    dist_rnorm_all   = vv2m( dist_rnorm_all ) 
    p1_Unorm_sum_all = vv2m( p1_Unorm_sum_all ) 
    p2_Unorm_sum_all = vv2m( p2_Unorm_sum_all ) 

    dist_norm_mean    = mean( dist_rnorm_all, dims = 1 )[:] 
    p1_Unorm_sum_mean = mean( p1_Unorm_sum_all, dims = 1 )[:] 
    p2_Unorm_sum_mean = mean( p2_Unorm_sum_all, dims = 1 )[:] 

    T = params.tof * params.k_tt_replan 

    # plot 
    fig = Figure( resolution = (600, 600) ) 

    title_string = string(params.strategy, " games = ", length(gameS), "\n", "mean cumsum norm of U vectors")
    ax1 = Axis( fig[1,1], xlabel = "time", title = title_string ) 
    tt = ( 0 : length(p1_Unorm_sum_mean)-1 ) * T / length(p1_Unorm_sum_mean)  
    p1_ax1 = lines!( ax1, tt, p1_Unorm_sum_mean, color = :blue ) 
    p2_ax1 = lines!( ax1, tt, p2_Unorm_sum_mean, color = :red ) 
    Legend( fig[1,2], [ p1_ax1, p2_ax1 ], ["p1", "p2"] ) 
    for ii in eachindex(gameS)
        lines!( ax1, tt, p1_Unorm_sum_all[ii,:][:], color = :blue, alpha = 0.1 ) 
        lines!( ax1, tt, p2_Unorm_sum_all[ii,:][:], color = :red, alpha = 0.1 ) 
    end 

    ax2 = Axis( fig[2,1], xlabel = "time (s)", title = "mean distance" ) 

    tt = ( 0 : length(dist_norm_mean)-1 ) * T / length(dist_norm_mean)  
    lines!( ax2, tt, dist_norm_mean, color = :green ) 
    for ii in eachindex(gameS)
        lines!( ax2, tt, dist_rnorm_all[ii,:][:], color = :green, alpha = 0.1 ) 
    end 

    return fig 
end 

export plot_MC_stats 

## ============================================ ##

function plot_game_stats( gameS, ii, params ) 

    game = gameS[ii] 

    r_norm = dist_norm( game, params ) 
    p1_U_norm, p2_U_norm = U_norm( game, params ) 
    p1_Unorm_sum = cumsum( p1_U_norm ) 
    p2_Unorm_sum = cumsum( p2_U_norm ) 

    fig = Figure( resolution = (600, 800) )

    title_string = string( params.strategy, " game = ", ii, "\n norm of U vectors" ) 
    ax1 = Axis( fig[1,1], xlabel = "time", title = title_string ) 
    p1_ax1 = lines!( ax1, 1 : length(p1_U_norm), p1_U_norm, color = :blue ) 
    p2_ax1 = lines!( ax1, 1 : length(p2_U_norm), p2_U_norm, color = :red ) 
    Legend( fig[1,2], [ p1_ax1, p2_ax1 ], ["p1", "p2"] ) 

    ax2 = Axis( fig[2,1], xlabel = "time", title = "cumsum of norm of U vectors" ) 
    lines!( ax2, 1 : length(p1_Unorm_sum), p1_Unorm_sum, color = :blue ) 
    lines!( ax2, 1 : length(p2_Unorm_sum), p2_Unorm_sum, color = :red ) 

    ax3 = Axis( fig[3,1], xlabel = "time", title = "distance" ) 
    lines!( ax3, 1 : length(r_norm), r_norm, color = :green )  

    return fig 
end 

export plot_game_stats 


