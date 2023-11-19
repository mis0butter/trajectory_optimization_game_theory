using GLMakie 

## ============================================ ## 

export plot_orbit 
function plot_orbit( 
    rv, 
    fig    = nothing, 
    labels = nothing 
    ) 

    if isnothing(fig) 
        # initialize figure 
        fig = Figure() 
        Axis3(fig[1, 1], 
            xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)", 
            title = "Transfer Solution") 
    end 

    text_offset = (0,10) 

    # plot solution 
    lines!( rv[:,1], rv[:,2], rv[:,3]; linewidth = 2, label = "lambert" ) 
    scatter!( rv[1,1], rv[1,2], rv[1,3]; marker = :utriangle, markersize = 15, color = :black ) 
    scatter!( rv[end,1], rv[end,2], rv[end,3]; marker = :xcross, markersize = 15, color = :black ) 
    if ~isnothing(labels) 
        text!( rv[1,1], rv[1,2], rv[1,3]; text = "start", color = :gray, offset = text_offset, align = (:center, :bottom) ) 
        text!( rv[end,1], rv[end,2], rv[end,3]; text = "end", color = :gray, offset = text_offset, align = (:center, :bottom) ) 
    end 

    return fig 
end 
    

