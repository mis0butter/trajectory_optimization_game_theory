using trajectory_optimization_game_theory
using LinearAlgebra 
using ForwardDiff 
using GLMakie 

## ============================================ ##

# define IC and target state 

mu = 398600.4415
r  = 6378.0
kep0_P = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_E = [r+450.0, 0.0, 10.6*pi/180, 0.0, 0.0, 50.0*pi/180]
t = (0.0, 1*orbitPeriod(kep0_E, mu)) 
t_P, x_P = propagate_2Body(kep2cart(kep0_P, mu), t, mu, 1.0)
t_E, x_E = propagate_2Body(kep2cart(kep0_E, mu), t, mu, 1.0)

x0_P_OG = x_P[1] 
x0_E_OG = x_E[1] 
xf_E_OG = x_E[end] 

## ============================================ ##

x0_P = copy(x0_P_OG)
x0_E = copy(x0_E_OG)  
xf_E = copy(xf_E_OG)  

x_P_dtsim_hist = [] ;     t_P_dtsim_hist = [] 
x_E_dtsim_hist = [] ;     t_E_dtsim_hist = [] 
x_P_tof_hist   = [] ;     t_P_tof_hist   = [] 
x_E_tof_hist   = [] ;     t_E_tof_hist   = [] 
dt_sim = 10 
k_sim  = 0 
for t_sim = dt_sim : dt_sim : 1000

    k_sim += 1 

    println(t_sim) 

    # t_sim = 10 
    # lambert 
    tof     = 100     # this can change / vary with function update 
    dm      = "retro" # replace with whatever is min 

    r1 = x0_P[1:3] 
    r2 = xf_E[1:3] 

    # solve transfer 
    v1 = lambertbattin(r1, r2, mu, dm, tof) 
    x0_P_lambert   = [r1; v1] 
    t_P_dtsim, x_P_dtsim = propagate_2Body( x0_P_lambert, dt_sim, mu, 1.0 )
    t_P_tof,   x_P_tof   = propagate_2Body( x0_P_lambert, tof, mu, 1.0 )
    xf_P = x_P_dtsim[end] 
    # for i in 1 : length(x_P_dtsim)
        x_P_dtsim = mapreduce( permutedims, vcat, x_P_dtsim ) 
        x_P_tof   = mapreduce( permutedims, vcat, x_P_tof ) 
        push!(x_P_dtsim_hist, x_P_dtsim)
        push!(x_P_tof_hist, x_P_tof)
    # end 

    # move evader 
    # x0_E += [0.0, 0.0, 0.0, 10/1000, 10/1000, 10/1000] 
    t_E_dtsim, x_E_dtsim = propagate_2Body(x0_E, dt_sim, mu, 1.0) 
    t_E_tof,   x_E_tof   = propagate_2Body(x0_E, tof, mu, 1.0) 
    xf_E = x_E_dtsim[end] 
    # for i in 1 : length(x_E_dtsim) 
        x_E_dtsim = mapreduce( permutedims, vcat, x_E_dtsim ) 
        x_E_tof   = mapreduce( permutedims, vcat, x_E_tof ) 
        push!(x_E_dtsim_hist, x_E_dtsim)
        push!(x_E_tof_hist, x_E_tof) 
    # end

    # update IC 
    x0_P = xf_P 
    x0_E = xf_E 

end 

# x_P_dtsim_hist = mapreduce( permutedims, vcat, x_P_dtsim_hist ) 
# x_E_dtsim_hist = mapreduce( permutedims, vcat, x_E_dtsim_hist ) 

# # plot 
# fig = plot_orbit( x_P_dtsim_hist )
# fig = plot_orbit( x_E_dtsim_hist, fig ) 
# fig = plot_scatter3d( x_P_dtsim_hist[1,1], x_P_dtsim_hist[1,2], x_P_dtsim_hist[1,3], fig ) 
# fig = plot_scatter3d( x_E_dtsim_hist[1,1], x_E_dtsim_hist[1,2], x_E_dtsim_hist[1,3], fig ) 

## ============================================ ##
# create animation 

# using DataFrames
# using GLMakie
using Plots 

a = Animation() 

for i = 2 : k_sim 

    x_P_dtsim = x_P_dtsim_hist[1]  
    x_E_dtsim = x_E_dtsim_hist[1] 
    for j = 2 : i 
        x_P_dtsim = [ x_P_dtsim ; x_P_dtsim_hist[j] ] 
        x_E_dtsim = [ x_E_dtsim ; x_E_dtsim_hist[j] ] 
        # push!( x_P_dtsim, x_P_dtsim_hist[j] ) 
        # push!( x_E_dtsim, x_E_dtsim_hist[j] ) 
    end 
    # x_P_dtsim = mapreduce( permutedims, vcat, x_P_dtsim ) 
    # x_E_dtsim = mapreduce( permutedims, vcat, x_E_dtsim ) 

    # x_P_dtsim = x_P_dtsim_hist[i] 
    # x_E_dtsim = x_E_dtsim_hist[i] 
    x_P_tof   = x_P_tof_hist[i] 
    x_E_tof   = x_E_tof_hist[i] 

    plt = plot3d(
        x_P_dtsim[:,1], x_P_dtsim[:,2], x_P_dtsim[:,3], 
        legend = false,
        color  = :blue, 
    ) 
    plot3d!(  
        x_P_tof[:,1], x_P_tof[:,2], x_P_tof[:,3], 
        linealpha = 0.5, 
        color  = :cyan, 
    ) 
    # scatter3d!(  
    #     [x_P_dtsim_hist[1,1]], [x_P_dtsim_hist[1,2]], [x_P_dtsim_hist[1,3]],
    #     # marker = 2, 
    #     markershape = :utriangle, 
    #     color  = :blue,
    # )
    # scatter3d!(  
    #     [x_P_dtsim_hist[i,1]], [x_P_dtsim_hist[i,2]], [x_P_dtsim_hist[i,3]],
    #     markershape = :xcross, 
    #     color  = :blue,  
    # )

    plot3d!(  
        x_E_dtsim[:,1], x_E_dtsim[:,2], x_E_dtsim[:,3], 
        color  = :red, 
    ) 
    plot3d!(  
        x_E_tof[:,1], x_E_tof[:,2], x_E_tof[:,3], 
        linealpha = 0.5, 
        color  = :orange, 
    )
    # scatter3d!(  
    #     [x_E_dtsim_hist[1,1]], [x_E_dtsim_hist[1,2]], [x_E_dtsim_hist[1,3]],
    #     markershape = :utriangle, 
    #     color  = :red,
    # )
    # scatter3d!(  
    #     [x_E_dtsim_hist[i,1]], [x_E_dtsim_hist[i,2]], [x_E_dtsim_hist[i,3]],
    #     markershape = :xcross, 
    #     color  = :red,  
    # )

    frame(a, plt)

end

g = gif(a, fps = 2.0)
display(g)  


