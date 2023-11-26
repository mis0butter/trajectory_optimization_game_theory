

#============================================================
SOLVE_TRANSFER:

Description: Solves the low-thrust transfer between an initial orbit state and a final orbit state using the Sims-Flanagan 
             method of discretization.


Inputs:
    xₒ   - Initial state vector constructed as [r;v]
    xfₒ  - Final state vector constructed as [r;v] 
    N    - Number of segments the trajectory is split up into
    μ    - Gravitational Parameter
    R    - Radius of the central body 
    revs - Number of revolutions between legs, initialized to 0


Output:
    Final Sims-Flanagan Trajectory object 
        • Look at Optimization.jl for the declaration of this struct
============================================================# 

function solve_transfer(
    x₀, 
    xfₒ, 
    N, 
    R, 
    μ;
    revs = 0
) 

    # Nondimensionalizing Inputs
    #   x̄₀ - Non-dimensionalized state vector
    #   DU - Relative distance
    #   TU - Relative time 
    x̄₀, DU, TU = nondimensionalize_x(x₀, μ, R) 

    # Other Initial Pieces
    #   xfₒ - final state at end of time period
    #   x̄f₀ - non-dimensionalized final state
    x̄f₀ = copy(xfₒ)
    x̄f₀[1:3] /= DU      
    x̄f₀[4:6] /= DU/TU
   
    # Number of Constraints
    mC = 3     # Terminal Constraints
    mI = N+1   # Number of Inequality Constraints

    # Initializing Performance Index
    ϕ = (x, λ, p) -> performance_index(x, λ, p, N, x̄₀, x̄f₀, mC)
    
    # Initializing Decision Variables
    #   ā - non-dimensionalized semi-major axis
    #   x_vec - vector of size 3N+1 containing Δτ and N velocity vectors for each segment of trajectory
    #   Δτ - Kepler's Universal Variable
    #   Δv_vec - N velocity vectors for each segment of trajectory 
    ā = x2oe(x̄₀)[1]
    x_vec = [sqrt(ā)*((revs + 1)*2π - π/8)/N, 1e-10*ones(3*N)...]
    Δτ, Δv_vec = unwrap(x_vec, N)

    # Initializing Pentalties
    p_vec = 10*ones(mC + mI)
    p_vec[mC+1:end] .= 0.0
    λ_vec = zeros(mC + mI)

    # Initializing Loop Variables
    γ = 2.0
    iter = 0
    isConverged = false
    while !isConverged && iter < 100
        iter += 1

        # Creating Function
        ϕ′  = x -> ϕ(x, λ_vec, p_vec)
        ∇ϕ′ = x -> ForwardDiff.gradient(ϕ′, x)

        # Main.@infiltrate 

        # Finding New Minimum
        x_new_vec = bfgs(x_vec, ϕ′, ∇ϕ′; tol = 1e-10, itermax=50)

        # Updating
        δx = norm(x_vec - x_new_vec)
        if δx < 1e-8
            isConverged = true
        end
        x_vec = copy(x_new_vec)

        # Evaluting Constraints
        Δτ, Δv_vec = unwrap(x_vec, N)
        h = equality_constraints(Δτ, Δv_vec, x̄₀, x̄f₀, N)
        ψ = inequality_constraints(Δτ, Δv_vec, x̄₀, x̄f₀, N)

        # Branching
        for i in eachindex(p_vec)

            # equality constraints 
            if i ≤ mC && abs(h[i]) ≥ 1e-4
                λ_vec[i] += p_vec[i]*abs(h[i])
                p_vec[i] *= γ

            # inequality constraints 
            elseif i > mC 

                if ψ[i - mC] > -λ_vec[i]/p_vec[i] 
                    λ_vec[i] = max(λ_vec[i] + p_vec[i]*ψ[i - mC], 0.0)
                    p_vec[i] *= ψ[i - mC] > 0.0 ? γ : 1
                end 

            end
        end

        if any(h .≥ 1e-3)
            isConverged = false
        end

        # Checking for Break Conditions
        if any(p_vec .≥ 1e12)
            @warn "Exitting due to high penalty values" p_vec maxlog = 1
            @info "Including Lagrange Multipliers" λ_vec maxlog = 1
            break
        end

        Δτ, Δv_vec = unwrap(x_vec, N)

    end
    Δv = [norm(Δv_vec[i, :]) for i in 1:N]
    x̄f, Δt̄ = prop_stateUV_Nseg(x̄₀, Δv_vec, Δτ, N)
    miss = equality_constraints(Δτ, Δv_vec, x̄₀, x̄f₀, N)
    
    xf = vcat(x̄f[1:3]*DU, x̄f[4:6]*DU/TU)
    Δt = Δt̄*TU

    return SimsFlanaganTrajectory(Δτ, Δv_vec*(DU/TU), Δv_vec, sum(Δv)*(DU/TU), 
        DU, TU, x₀, xf, Δt, miss)
end


#============================================================
UNWRAP:

Description: Helper function for solve_transfer that breaks up a specific type of vector into two different components

Inputs:
    1. x_vec - the vector to be broken down
        •   has the form of [Δτ, Δv₁, Δv₂, ..., Δvₙ]
    2. N - Number of segments of the trajectory

Outputs:
    1. Δτ - Kepler's Universal Variable
    2. Δv_vec - vector of size (N,3) containing velocity vectors for each segment of the trajectory

============================================================#

function unwrap(
    x_vec::AbstractVector, 
    N::Int)

    Δτ = x_vec[1]
    Δv_vec = reshape(x_vec[2:end], N, 3)

    return Δτ, Δv_vec
end
