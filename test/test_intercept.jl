using trajectory_optimization_game_theory
using LinearAlgebra 
using ForwardDiff 
using GLMakie 

## ============================================ ##

# define IC and target state 

mu = 398600.4415
r  = 6378.0
kep0_P = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_E = [r+400.0, 0.0, 10.6*pi/180, 0.0, 0.0, 50.0*pi/180]
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

x_P_hist = [] ;     t_P_hist = [] 
x_E_hist = [] ;     t_E_hist = [] 
dt_sim = 10 
for t_sim = dt_sim : dt_sim : 1000 

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
    t_P_lambert, x_P_lambert = propagate_2Body( x0_P_lambert, dt_sim, mu, 1.0 )
    xf_P = x_P_lambert[end] 
    for i in 1 : length(x_P_lambert)
        push!(x_P_hist, x_P_lambert[i])
    end 

    # move evader 
    # x0_E += [0.0, 0.0, 0.0, 10/1000, 10/1000, 10/1000] 
    t_E, x_E = propagate_2Body(x0_E, dt_sim, mu, 1.0) 
    xf_E = x_E[end] 
    for i in 1 : length(x_E) 
        push!(x_E_hist, x_E[i])
    end

    # update IC 
    x0_P = xf_P 
    x0_E = xf_E 

end 

x_P_hist = mapreduce( permutedims, vcat, x_P_hist ) 
x_E_hist = mapreduce( permutedims, vcat, x_E_hist ) 

# plot 
fig = plot_orbit( x_P_hist )
fig = plot_orbit( x_E_hist, fig ) 
fig = plot_scatter3d( x_P_hist[1,1], x_P_hist[1,2], x_P_hist[1,3], fig ) 
fig = plot_scatter3d( x_E_hist[1,1], x_E_hist[1,2], x_E_hist[1,3], fig ) 

## ============================================ ##
# create animation 

# using DataFrames
# using GLMakie
using Plots 

a = Animation() 

for i = 2 : size(x_P_hist, 1)

    plt = plot3d(
        x_P_hist[1:i,1], x_P_hist[1:i,2], x_P_hist[1:i,3], 
        legend = false,
        color  = :blue, 
    ) 
    scatter3d!(  
        [x_P_hist[1,1]], [x_P_hist[1,2]], [x_P_hist[1,3]],
        # marker = 2, 
        markershape = :utriangle, 
        color  = :blue,
    )
    scatter3d!(  
        [x_P_hist[i,1]], [x_P_hist[i,2]], [x_P_hist[i,3]],
        markershape = :xcross, 
        color  = :blue,  
    )

    plot3d!(  
        x_E_hist[1:i,1], x_E_hist[1:i,2], x_E_hist[1:i,3], 
        color  = :red, 
    ) 
    scatter3d!(  
        [x_E_hist[1,1]], [x_E_hist[1,2]], [x_E_hist[1,3]],
        markershape = :utriangle, 
        color  = :red,
    )
    scatter3d!(  
        [x_E_hist[i,1]], [x_E_hist[i,2]], [x_E_hist[i,3]],
        markershape = :xcross, 
        color  = :red,  
    )

    frame(a, plt)

end

g = gif(a, fps = 10.0)
display(g)  


