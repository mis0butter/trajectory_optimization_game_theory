using GLMakie 

## ============================================ ## 

"Plot an orbit using GLMakie"
function plot_orbit( 
    rv,                 # [N,3] matrix of state vectors 
    fig    = nothing,   # figure handle 
    labels = false      # boolean for labeling start and end points 
) 

    if isnothing(fig) 
        fig = Figure() 
        Axis3(fig[1, 1], 
            xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)", 
            title = "Transfer Solution") 
    end 

    text_offset = (0,10) 

    # plot orbit 
    lines!( rv[:,1], rv[:,2], rv[:,3]; linewidth = 2, label = "lambert" ) 
    scatter!( rv[1,1], rv[1,2], rv[1,3]; marker = :utriangle, markersize = 15, color = :black ) 
    scatter!( rv[end,1], rv[end,2], rv[end,3]; marker = :xcross, markersize = 15, color = :black ) 

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

"Plot a surface with a colorbar using GLMakie"
function plot_surface( 
    x,              # [N,1] grid of points 
    y,              # [N,1] grid of points 
    z,              # [N,N] grid of points evaluated at x and y 
    fig = nothing   # figure handle 
) 

    fignothing = false 
    if isnothing(fig) 
        fignothing = true 
        fig = Figure() 
        Axis3( fig[1,1] ) 
    end 

    if isequal(length(z), 1)
        GLMakie.scatter!( x, y, z, marker = :utriangle, markersize = 20, color = :red ) 
    else 
        hm = GLMakie.surface!( x, y, z, colormap = :diverging_tritanopic_cwr_75_98_c20_n256 ) 
    end 

    if fignothing 
        Colorbar( fig[1,2], hm, height = Relative(0.5) )
    end 

    return fig 
end 

export plot_surface 

## ============================================ ##

"Plot a contour with a colorbar using GLMakie"
function plot_contour3d( 
    x,              # [N,1] grid of points 
    y,              # [N,1] grid of points 
    z,              # [N,N] grid of points evaluated at x and y 
    fig = nothing,  # figure handle 
    levels = 20,    # number of contour levels 
) 

    fignothing = false 
    if isnothing(fig) 
        fignothing = true 
        fig = Figure() 
        Axis3(fig[1,1]) 
    end 

    hm  = GLMakie.contour3d!(x, y, z, levels = levels) 

    if fignothing 
        clim = ( minimum(z), maximum(z) ) 
        Colorbar( fig[1, 2], limits = clim, height = Relative(0.5) )
    end 

    return fig 
end 

export plot_contour3d 


