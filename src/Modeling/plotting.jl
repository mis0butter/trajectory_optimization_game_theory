using GLMakie 

#============================================================
PLOT_SOLUTION!:
 
Description: Plots the trajectory between x₀ and xf₀

Inputs:
    x₀      - Initial state vector
    xf₀     - Final state vector
    Δτ      - Kepler's Universal Variable
    Δv_vec  - Matrix of size (N,3) where each row is a velocity vector at a segment of the trajectory
    μ       - Gravitational Parameter
    label   - Label of graph, initialized to nothing
    color   - Color of the line of the graph, initialized to nothing (results in blue line)

Outputs:
    fig     - Figure object, displays graph
============================================================# 

# using Infiltrator 

function plot_solution!(
    x₀ ,
    xf₀,
    Δτ,
    Δv_vec, 
    R, 
    μ  = 1.0,
    label = nothing, 
    color = nothing, 
    fig   = nothing
    ) 
    
    # Non-DimensionalizingS
    x₀, DU, TU = nondimensionalize_x(x₀, μ, R)
    xf₀ = copy(xf₀)
    xf₀[1:3]  /= DU
    xf₀[4:6]  /= DU/TU
    Δv_vec    /= DU/TU

    # Getting Required Trajectory States
    N = size(Δv̄_vec, 1)
    Xtraj, Δt = prop_stateUV_Nseg_range(x̄₀, Δv̄_vec, Δτ, 1:N)

    # println( "integrated Xtraj" ) 
    # @infiltrate 

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

    # println( "integrated Xtraj2" ) 
    # @infiltrate 

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

    # println( "Redimmensionalizing" ) 
    # @infiltrate 

    # Initializing Figure
    if isnothing(fig)  
        fig = Figure()
        Δv = sum([norm(Δv̄_vec[i, :]) for i in 1:N])
        error = norm(Xtraj[end, 1:3] - Xf[end, 1:3])
        vinf  = norm(Xtraj[end, 4:6] - Xf[end, 4:6])
         Axis3(fig[1, 1], 
             xlabel = "X (km)", ylabel = "Y (km)", zlabel = "Z (km)", 
             title = "Total Δv of Transfer: $(round(Δv; sigdigits=6)*1000) m/s\nFinal Error: $(round(error; sigdigits=6)) km\nVinf: $(round(vinf; sigdigits=6)) km/s")
    end 

    # Plotting Initial Condition
    # Plots line between the initial state and the propagated initial state
     lines!(Xi[:, 1], Xi[:, 2], Xi[:, 3]; 
         linestyle =  :dash, color = :black)

    # println( "plotted IC" ) 
    # @infiltrate 
     
    # Plotting Final Condition
    # Plots line between the final state and the propagated final state
     lines!(Xf[:, 1], Xf[:, 2], Xf[:, 3]; 
        linestyle =  :dot, color = :black)

    # println( "plotted final condition" ) 
    # @infiltrate     

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

    Auto() 

    # Reshowing Figure
    return fig
end

#============================================================
PLOTTRAJ:

Description: Plots a trajectory, approximating between states using Lambert

Inputs:
    1. Δv_vec - Matrix of size (N, 3) where each row is the velocity vector at a segment of the trajectory
    2. Δτ - Kepler's Universal Variable
    3. x0 - Initial state
    4. xf0 - Final state
    5. μ - Gravitational parameter

Outputs:
    1. Plot of the trajectory
============================================================#
function plottraj(Δv_vec, Δτ, x0, xf0, μ, R)
    x̄0, DU, TU = UpdatedAdversarialTourDesign.nondimensionalize_x(x0, μ, R)
    # Getting Trajectory States
    X = getstates(Δv_vec, Δτ, x0, μ)

    # Initializing Figure
    fig = Figure()
    ax = Axis3(fig[1, 1])

    # Propagating Initial State
    a0, _, _, _, _, _ = x2oe(x0, μ)
    T = 2π*sqrt(a0^3 / μ)
    N = 100
    X0 = zeros(6, N)
    δt = LinRange(0, T, N)
    for i = 1:N
        X0[:, i] = propagate_x(xf0, δt[i], μ)
    end
    lines!(X0[1, :]/DU, X0[2, :]/DU, X0[3,:]/DU;
       linestyle = :dash)

    # Initial Guess Lambert
    xf = propagate_x(xf0, T, μ)
    l = lambert(x0, xf, T, μ)
    @show l.Δvᵢ
    Xl = zeros(6, N)
    δt = LinRange(0, T, N)
    for i = 1:N
        Xl[:, i] = propagate_x(l.xᵢ, δt[i], μ)
    end

    # Plotting Guess
    lines!(Xl[1, :]/DU, Xl[2, :]/DU, Xl[3,:]/DU;
        linestyle = :dash)

    # Plotting Trajectory
    @info "Final State" X[1:3, end]/DU
    lines!(X[1, :]/DU, X[2, :]/DU, X[3,:]/DU;
        linewidth=2)

    # Plotting dvs
    arrows!(X[1, 1:end-1]/DU, X[2, 1:end-1]/DU, X[3, 1:end-1]/DU, Δv_vec[:, 1], Δv_vec[:, 2], Δv_vec[:,3];
    linewidth = .1, arrowsize = .1)

    fig
end

#============================================================
GETSTATES:

Decription: Helper function for plottraj that calculates the state at each segment of the trajectory

Inputs: 
    1. Δv_vec - Matrix of size (N, 3) where each row is the velocity vector at a segment of the trajectory
    2. Δτ - Kepler's Universal Variable
    3. x0 - Initial state
    4. μ - Gravitational parameter
============================================================#

function getstates(Δv_vec, Δτ, x0, μ, R)
    # Non-Dimensionalizing State
    @show x0, DU, TU = nondimensionalize_x(x0, μ, R)
    # @show Δv⃗[1]
    Δv_vec /= (DU/TU)
    # @show Δv⃗[1]

    # Number of Δv's
    N = size(Δv_vec, 1)

    # Preallocating States
    X = zeros(6, N+1)
    X[:, 1] = x0

    # Propagating to Final State
    x  = x0
    Δt = 0.0
    Δx = zeros(6)
    for i = 1:N
        # Applying Δv at **beginning** of each node
        # @info "test" x[4:6] Δv⃗[i, :] maxlog = 1
        
        Δx[4:6] = Δv_vec[i, :]
        # @show Δx
        x += Δx

        # Propagating Forward Using Universial Variable
        x, δt = propKepUV(x, Δτ*10)
        Δt += δt
        X[:, 1+i] = x
    end

    X[1:3, :] *= DU
    X[4:6, :] *= DU/TU

    return X
end