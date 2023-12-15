using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 

## ============================================ ##
# init params 

mu = 398600.4415
r  = 6378.0
kep0_P = [ r+400.0, 0.1, -20*pi/180, 10.0*pi/180, 20.0*pi/180, 30.0*pi/180 ]
rv_0_P = rv_0_P_OG = kep2cart(kep0_P, mu) 
kep0_E = [ r+450.0, 0.2, 10.6*pi/180, 40.0*pi/180, 0.0, 40.0*pi/180 ]
rv_0_E = rv_0_E_OG = kep2cart(kep0_E, mu) 

# tof for pursuer to catch up to evader 
tof = 2000 

t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0) 
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0) 
rv_P = vv2m(rv_P) 
rv_E = vv2m(rv_E) 

# plot 
fig = plot_axes3d( )
fig = plot_orbit( rv_P, fig ) 
fig = plot_orbit( rv_E, fig ) 

## ============================================ ##
# break up into N segments, see what happens 

# define init and target vectors for pursuer 
rv_f = rv_E[end,:] 
rv_0 = rv_0_P 

N = 20 
tof_N = tof / N 
Δv_sol = min_Δv( rv_0, rv_f, tof, N, mu ) 
Δv_sol2 = min_Δv_dist( rv_0, rv_f, tof, N, mu ) 

t, rv_kepler = prop_kepler_tof_Nseg( rv_0, Δv_sol, N, tof / N, mu ) 

## ============================================ ##
# MPC 

# init 
rv_0_E = copy( rv_0_E_OG ) 
rv_0_P = copy( rv_0_P_OG ) 

rv_P_dtsim_hist = [] ;  t_P_dtsim_hist = [] 
rv_E_dtsim_hist = [] ;  t_E_dtsim_hist = [] 
rv_P_tof_hist = [] ;    t_P_tof_hist = [] 
rv_E_tof_hist = [] ;    t_E_tof_hist = [] 

N_prop = 2 
dt_sim = tof_N * N_prop 
k_sim  = 0 
for t_sim = dt_sim : dt_sim : dt_sim * 5 

    k_sim += 1 
    println( "k_sim = ", k_sim ) 

    # propagate evader based on current information for tof 
    t_E_tof, rv_E_tof = propagate_2Body( rv_0_E, tof, mu, 1.0 ) 
    push!( rv_E_tof_hist, rv_E_tof )  
    push!( t_E_tof_hist, t_E_tof ) 

    rv_f_E   = rv_E_tof[end] 
    rv_E_tof = vv2m( rv_E_tof ) 

    # solve for pursuer Δv for tof 
    Δv_P = min_Δv( rv_0_P, rv_f_E, tof, N, mu ) 

    # compute entire predicted trajectory for pursuer  
    # t_P_tof, rv_P_tof = prop_kepler_tof_Nseg( rv_0_P, Δv_P, N, tof_N, mu ) 
    t_P_tof, rv_P_tof = prop_2Body_tof_Nseg( rv_0_P, Δv_P, N, tof_N, mu ) 
    push!( rv_P_tof_hist, rv_P_tof ) 
    push!( t_P_tof_hist, t_P_tof ) 

    # propagate pursuer for N = 10 segments 
    # t_P_dtsim, rv_P_dtsim = prop_kepler_tof_Nseg( rv_0_P, Δv_P, N_prop, tof_N, mu ) 
    t_P_dtsim, rv_P_dtsim = prop_2Body_tof_Nseg( rv_0_P, Δv_P, N_prop, tof_N, mu ) 
    push!( rv_P_dtsim_hist, rv_P_dtsim ) 
    push!( t_P_dtsim_hist, t_P_dtsim ) 

    # set up for next iter 
    rv_0_P = rv_P_dtsim[end,:] 

    # apply maneuver to evader 
    Δi = 10.0*pi/180 
    Δv_E = computeInclinationChange(rv_0_E, Δi, mu)
    rv_0_E += [0.0, 0.0, 0.0, Δv_E[1], Δv_E[2], Δv_E[3]]

    # propagate evader for dt_sim 
    t_E_dtsim, rv_E_dtsim = propagate_2Body( rv_0_E, t_sim, mu, 1.0 ) 
    rv_E_dtsim = vv2m( rv_E_dtsim ) 
    push!( rv_E_dtsim_hist, rv_E_dtsim )  
    push!( t_E_dtsim_hist, t_E_dtsim ) 

    # set up for next iter 
    rv_0_E = rv_E_dtsim[end,:] 


end 


## ============================================ ##
# create animation 

using Plots 
# using Formatting 

plot_font = "Computer Modern"
default(
    fontfamily = plot_font,
    linewidth  = 2, 
    framestyle = :box, 
    label      = nothing, 
    grid       = false, 
    markerstrokewidth = 0, 
)
# scalefontsizes(1.3)

# using DataFrames
# using GLMakie 

a = Animation() 

for i = 2 : k_sim 

    # get total history of x dtsim as matrix 
    rv_P_dtsim = rv_P_dtsim_hist[1]  
    rv_E_dtsim = rv_E_dtsim_hist[1] 
    for j = 2 : i 
        rv_P_dtsim = [ rv_P_dtsim ; rv_P_dtsim_hist[j] ] 
        rv_E_dtsim = [ rv_E_dtsim ; rv_E_dtsim_hist[j] ] 
    end 

    # get current iterate of x tof as matrix 
    rv_P_tof   = rv_P_tof_hist[i] 
    rv_E_tof   = rv_E_tof_hist[i] 

    # set up plot 
    lims = 1.1 .* (-r, r)
    plt = plot3d(
        1,
        # xlim = lims,
        # ylim = lims,
        # zlim = lims,
        title = "Pursuit Evasion Game",
        legend = false,
        xlabel = "x (km)", 
        guidefontsize = 10,
        tickfontsize = 8,
        ylabel = "y (km)",
        zlabel = "z (km)",
        titlefont=font(16,"Computer Modern"), 
        # camera = (-30,35,30), 
    )
    
    # Plots.surface!( 
    #     sphere(r, zeros(3)), 
    #     alpha = 0.1 
    # )

    # Pursuer!!!     
    plot3d!(
        rv_P_dtsim[:,1], rv_P_dtsim[:,2], rv_P_dtsim[:,3], 
        legend = false,
        color  = :blue, 
    ) 
    plot3d!(  
        rv_P_tof[:,1], rv_P_tof[:,2], rv_P_tof[:,3], 
        linealpha = 0.5, 
        color  = :cyan, 
    ) 
    scatter3d!(  
        [rv_P_dtsim[1,1]], [rv_P_dtsim[1,2]], [rv_P_dtsim[1,3]],
        # marker = 2, 
        markershape = :xcross, 
        color  = :blue,
    )
    scatter3d!(  
        [rv_P_dtsim[end,1]], [rv_P_dtsim[end,2]], [rv_P_dtsim[end,3]],
        markershape = :utriangle, 
        color  = :blue,  
    )

    # Evader!!! 
    plot3d!(  
        rv_E_dtsim[:,1], rv_E_dtsim[:,2], rv_E_dtsim[:,3], 
        color  = :red, 
    ) 
    plot3d!(  
        rv_E_tof[:,1], rv_E_tof[:,2], rv_E_tof[:,3], 
        linealpha = 0.5, 
        color  = :orange, 
    )
    scatter3d!(  
        [rv_E_dtsim[1,1]], [rv_E_dtsim[1,2]], [rv_E_dtsim[1,3]],
        markershape = :xcross, 
        color  = :red,
    )
    scatter3d!(  
        [rv_E_dtsim[end,1]], [rv_E_dtsim[end,2]], [rv_E_dtsim[end,3]],
        markershape = :utriangle, 
        color  = :red,  
    )

    frame(a, plt)

    filename_string = string( "test/outputs/test_intercept_", i, ".png" )
    savefig(filename_string)

end

g = gif(a, fps = 2.0)
display(g)  




