#============================================================
BFGS:

Description: Unconstrained gradient descent solver using the Broyden-Fletcher-Goldfarb-Shanno algorithm.

Inputs:
    1. x₀ - Initial state vector
    2. f - function of x
        •   Example: solve_transfer() uses a performance index function
    3. g - Derivative of f with respect to x 
    4. tol - Convergence tolerance, initialized to 1e-8
    5. itermax - Maximum amount of iterations
    6. params - additional parameters for f and g

Outputs:
    1. x - Final state vector
============================================================#

function bfgs(x₀, f, g; tol = 1e-8, itermax = 100, params = (;))

    # Initializing
    δx = Inf
    n = length(x₀)
    gk = zeros(n)
    Qk = zeros(n, n)
    xk = zeros(n)

    # Anonymizing Function
    if ~isempty(params)
        f = x -> f(x, params)
        g = x -> g(x, params)
    end

    # Iterating
    x = copy(x₀)
    iter = 0
    miter = 0
    exit = 0
    while abs(δx) > tol
        iter += 1

        # Finding Update Direction
        ŝ = update_gradients!(x, g, gk, Qk, xk)

        # Linesearching in that Direction
        bounds, lsiter = linesearch(x, ŝ, f)

        # Locating Minimum
        xnew, griter = goldenRatioMin(bounds, f)

        # Updating
        δx = norm(x - xnew)
        miter += lsiter + griter

        if iter == itermax 
            x = copy(xnew)
            exit = 1
            break
        end
        x = copy(xnew)
    end 

    # Outputting
    possible_info = [
        "Solution converged successfully!",
        "Solution did not converge: Ran out of iterations"
    ]

    if exit == 1
        @info possible_info[2]
    else
        @info possible_info[1]
    end
        
    return x
end

#============================================================
GOLDENRATIOMIN:

Description: Uses the golden ratio to find the minimum x value between two bounds

Inputs:
    1. x - bounds for the state vectors
    2. f - function of x 
    3. tol - Convergence tolerance
    4. itermax - Maximum number of iterations

Outputs:
    1. xf - Minimum x value
    2. iter - Number of iterations
============================================================#

function goldenRatioMin(x, f, tol = 1e-8, itermax = 50)

    # Initializing
    ϕ  = (3 - √5)/2
    xL = x[1]            ; fL = f(xL)
    xR = x[3]            ; fR = f(xR)
    x1 = xL + ϕ*(xR - xL); f1 = f(x1)
    x2 = xR - ϕ*(xR - xL); f2 = f(x2)
    δx = Inf
    iter = 0

    # Iterating
    while abs(δx) > tol
        iter += 1

        # Branching
        if f1 > f2
            # Adjusting Bounded Points
            xL = x1; fL = f1;
            x1 = x2; f1 = f2;
            
            # Re-Evaluating New Points
            x2 = xR - ϕ*(xR - xL);
            f2 = f(x2);
        else
            # Adjusting Bounded Points
            xR = x2; fR = f2;
            x2 = x1; f2 = f1;
            
            # Re-Evaluating New Points
            x1 = xL + ϕ*(xR - xL);
            f1 = f(x1);
        end

        # Finding Bounds
        δx = norm(xR - xL)

        if iter == itermax; 
            @warn "GoldRatioMin iteration failure"
            break
        end
    end

    # Finding Best Point
    X = [xL, x1, x2, xR]
    Z = [fL, f1, f2, fR]
    idx = findmin(Z)[2]

    xf = X[idx]

    return xf, iter
end

#============================================================
LINESEARCH:

Description: Line search algorithm

Inputs:
    1. x₀ - Initial state vector
    2. ŝ - gradient vector of f
    3. f - function of x
    4. initial_step - initial step size, initialized at 0.00001
    5. updatestep - factor saying how often to update the step size, initialized to 4
        •   i.e how many iterations before the step size is updated
    6. itermax - Maximum iterations
    7. stepmod - Factor by which the step size is updated
============================================================#

function linesearch(
    x₀::AbstractVector{T}, 
    ŝ::AbstractVector{T}, 
    f::Function; 
    initial_step::T = 0.00001, 
    updatestep::Int = 4, 
    itermax::Int = 50,
    stepmod::T = 2.0) where T<:AbstractFloat

    # Initializing
    z = [f(x₀)]
    t = initial_step
    x = [x₀]

    # Checking Initial Step Size
    if f(x₀ + t*ŝ) > z[1]
        @warn "First step increasing: StepSizeTooLarge" initial_step
        return [x₀, x₀ + 0.5*t*ŝ, x₀ + t*ŝ], 0
    end

    # Bracketting Minimum
    iter = 0
    while length(z) < 3 || z[end] < z[end-1]
        # Updating Iteration and Step Size
        iter += 1
        if mod(iter, updatestep) == 0
            t *= stepmod
        end

        # Pathing In Direction
        push!(x, x[end] + t*ŝ)
        push!(z, f(x[end]))

        if iter == itermax; 
            @warn "LineSearch iteration failure"
            break
        end
    end

    return x[end-2:end], iter
end

#============================================================
UPDATE_GRADIENTS:

Description: Function that solves for the gradients of a function at an x value

Inputs:
    1. x - Initial state vector
    2. g - Derivative function of f, a function of x 
    3. glast - initial guess of the derivative
    4. Qlast - initial guess of the Hessian
    5. xlast - values of x used to modify the Hessian

Outputs:
    1. ŝ - gradient vector
============================================================#

function update_gradients!(
    x::AbstractVector{T}, 
    g::Function, 
    glast::AbstractVector{T},
    Qlast::AbstractMatrix{T},
    xlast::AbstractVector{T}) where T<:AbstractFloat


    n = length(x)
    if iszero(glast)
        gk1 = g(x)
        Qk1 = Matrix(I, n, n)

    else
        # Importing
        gk1 = g(x)
        Qk = copy(Qlast)

        # Constants
        p = x - xlast
        y = gk1 - glast
        σ = p'*y
        τ = y'*Qk*y
        A = Qk*y*p'

        # Calcuating Approximated Hessian Inverse
        δQ = ((σ + τ)/σ^2)*(p*p') - (1/σ)*(A + A')
        Qk1 = Qk + δQ
    end

    # Updating Inputs
    xlast[:] = x
    glast[:] = gk1
    Qlast[:] = Qk1

    # Outputting Search Direction
    ŝ = -Qk1*gk1 / norm(Qk1*gk1)
    return ŝ
end