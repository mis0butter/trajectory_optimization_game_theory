using trajectory_optimization_game_theory
using LinearAlgebra, StaticArrays
using Infiltrator

# Infiltrator.clear_disabled!()
# Infiltrator.toggle_async_check(false) 

## ============================================ ##

mu = 398600.4415
r = 6378.0
kep0_p1 = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_p2 = [r+450.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]
t = (0.0, orbitPeriod(kep0_p2, mu))
prop1 = propagate_2Body(kep2cart(kep0_p1, mu), t, mu)
prop2 = propagate_2Body(kep2cart(kep0_p2, mu), t, mu)

x₀_P = prop1.u[1] 
x₀_E = prop2.u[1] 
xf_E = prop1.u[end] 

xₒ  = x₀_P 
xfₒ = xf_E 

sf  = solve_transfer(xₒ, 20, xfₒ, 0.0, mu, r)

## ============================================ ##

fig = plot_solution!(xₒ, sf.xf, sf.Δτ, sf.Δv⃗, mu) 

## ============================================ ##

x₀      = xₒ 
xf₀     = sf.xf 
Δτ      = sf.Δτ 
Δv_vec  = sf.Δv⃗ 
μ       = mu 
label   = nothing 
color   = nothing 
fig     = nothing 

## ============================================ ##



#============================================================
PLOT_SOLUTION!:

Description: Plots the trajectory between x₀ and xf₀

Inputs:
    1. x₀ - Initial state vector
    2. xf₀ - Final state vector
    3. Δτ - Kepler's Universal Variable
    4. Δv_vec - Matrix of size (N,3) where each row is a velocity vector at a segment of the trajectory
    5. μ - Gravitational Parameter
    6. label - Label of graph, initialized to nothing
    7. color - Color of the line of the graph, initialized to nothing (results in blue line)

Outputs:
    1. fig - Figure object, displays graph
============================================================# 

using GLMakie 

# using Infiltrator 

# function plot_solution!(
#     x₀ ::AbstractVector{T},
#     xf₀::AbstractVector{T},
#     Δτ::T,
#     Δv_vec::AbstractMatrix{T},
#     μ::T  = 1.0,
#     label = nothing, 
#     color = nothing, 
#     fig   = nothing) where T<:AbstractFloat 

    # Non-DimensionalizingS
    x̄₀, DU, TU = trajectory_optimization_game_theory.nondimensionalize_x(x₀, μ, r)
    x̄f₀    = copy(xf₀)
    Δv̄_vec = copy(Δv_vec)
    x̄f₀[1:3] /= DU
    x̄f₀[4:6] /= DU/TU 
    Δv̄_vec   /= DU/TU 

    # Getting Required Trajectory States
    N = size(Δv̄_vec, 1)
    Xtraj, Δt = trajectory_optimization_game_theory.prop_stateUV_Nseg_range(x̄₀, Δv̄_vec, Δτ, 1:N)

    println( "integrated Xtraj" ) 
    @infiltrate 

    # Integrating between segments
    m = 20
    Xtraj2 = Vector{Float64}[]
    for i = 1:N
        X0 = Xtraj[i, :]
        dτ = LinRange(0, Δτ, m)
        for j = 1:m
            push!(Xtraj2, propKepUV(X0, dτ[j])[1])
        end
    end
    Xtraj2 = hcat(Xtraj2...)'

    println( "integrated Xtraj2" ) 
    @infiltrate 

    # Propagating Orbit of initial body
    tspan = LinRange(0, Δt, N)
    Xi = zeros(N, 6)
    for i = 1:N
        Xi[i, :] = propagate_x(x̄₀, tspan[i],μ)
    end

    # Propagating Orbit of target body
    tspan = LinRange(0, Δt, N)
    Xf = zeros(N, 6)
    for i = 1:N
        Xf[i, :] = propagate_x(x̄f₀, tspan[i],μ)
    end

    # Redimmensionalizing
    Xtraj[:, 1:3] *= DU
    Xtraj[:, 4:6] *= DU/TU
    Xtraj2[:, 1:3] *= DU
    Xtraj2[:, 4:6] *= DU/TU
    Xi[:, 1:3] *= DU
    Xi[:, 4:6] *= DU/TU
    Xf[:, 1:3] *= DU
    Xf[:, 4:6] *= DU/TU
    Δv̄_vec         *= DU/TU

    println( "Redimmensionalizing" ) 
    @infiltrate 

    # Initializing Figure
    if isnothing(fig)  
        fig = Figure()
        Δv = sum([norm(Δv̄_vec[i, :]) for i in 1:N])
        err = norm(Xtraj[end, 1:3] - Xf[end, 1:3])
        vinf = norm(Xtraj[end, 4:6] - Xf[end, 4:6])
         Axis3(fig[1, 1], 
             xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)", 
             title = "Total Δv of Transfer: $(round(Δv; sigdigits=6)*1000) m/s\nFinal Error: $(round(err; sigdigits=6)) km\nVinf: $(round(vinf; sigdigits=6)) km/s")
    end 

    # Plotting Initial Condition
    # Plots line between the initial state and the propagated initial state
     lines!(Xi[:, 1], Xi[:, 2], Xi[:, 3]; 
         linestyle =  :dash, color = :black)

    println( "plotted IC" ) 
    @infiltrate 
     
    # Plotting Final Condition
    # Plots line between the final state and the propagated final state
     lines!(Xf[:, 1], Xf[:, 2], Xf[:, 3]; 
        linestyle =  :dot, color = :black)

    println( "plotted final condition" ) 
    @infiltrate     

    # Plotting Trajectory
    if isnothing(color)
        lines!(Xtraj2[:, 1], Xtraj2[:, 2], Xtraj2[:, 3]; 
            linewidth = 2)
    else
        if isnothing(label)
            lines!(Xtraj2[:, 1], Xtraj2[:, 2]. Xtraj2[:, 3]; color = Cycled(color),
                linewidth = 2)
        else
            lines!(Xtraj2[:, 1], Xtraj2[:, 2], Xtraj2[:, 3]; color = Cycled(color),
                linewidth = 2, label = label)
        end
    end
    scatter!(Xtraj[1, 1], Xtraj[1, 2], Xtraj[1,3]; 
        marker = :circle, markersize = 10, color = :black)
    scatter!(Xtraj[end, 1], Xtraj[end, 2], Xtraj[end,3]; 
        marker = :rect, markersize = 10, color = :black)

    # # Plotting DVs
    # Δvmax = maximum([norm(Δv_vec[i, :]) for i in 1:N])
    #  arrows!(Xtraj[1:end-1, 1], Xtraj[1:end-1, 2], Xtraj[1:end-1, 3], 
    #      Δv_vec[:, 1], 
    #      Δv_vec[:, 2], 
    #      Δv_vec[:, 3];
    #      lengthscale = 1e4/Δvmax, 
    #      linewidth = 300, 
    #      arrowsize = 200  
    #      )

    scatter!(0,0,0;marker = :circle, markersize = 15, color = :black)
    text!(0,0,0; text = "Earth?", color = :gray, offset = (-25,-25), align = (:center, :bottom))

    # Auto() 
    save("plot3d.png", fig)

    # Reshowing Figure
#     return fig
# end

## ============================================ ##





