""" Write a standard-form constrained static game as a MCP
i.e., if player i's problem is of the form
             min_{xᵢ} fᵢ(x, θ)
             s.t.     gᵢ(x, θ) = 0
                      hᵢ(x, θ) ≥ 0

and players share constraints
                      g̃(x, θ) = 0
                      h̃(x, θ) ≥ 0

Here, xᵢ is Pi's decision variable, x is the concatenation of all players'
variables, and θ is a vector of parameters. Express it in the following form:
             find   z
             s.t.   F(z, θ) ⟂ z̲ ≤ z ≤ z̅

where we interpret z = (x, λ, μ, λ̃, μ̃), with λ and μ the Lagrange multipliers
for the constraints g and h (with tildes as appropriate), respectively. The
expression F(z) ⟂ z̲ ≤ z ≤ z̅ should be read as the following three statements:
            - if z = z̲, then F(z, θ) ≥ 0
            - if z̲ < z < z̅, then F(z, θ) = 0
            - if z = z̅, then F(z, θ) ≤ 0

For more details, please consult the documentation for the package
`Complementarity.jl`, which may be found here:
https://github.com/chkwon/Complementarity.jl/tree/master
"""

"Generic description of a constrained parametric game problem."
struct ParametricGame{T1,T2,T3,T4,T5,T6,T7}
    "Objective functions for all players"
    objectives::T1
    "Equality constraints for all players"
    equality_constraints::T2
    "Inequality constraints for all players"
    inequality_constraints::T3
    "Shared equality constraint"
    shared_equality_constraint::T4
    "Shared inequality constraint"
    shared_inequality_constraint::T5

    "Dimension of parameter vector"
    parameter_dimension::T6
    "Dimension of primal variables for all players"
    primal_dimensions::T7
    "Dimension of equality constraints for all players"
    equality_dimensions::T7
    "Dimension of inequality constraints for all players"
    inequality_dimensions::T7
    "Dimension of shared equality constraint"
    shared_equality_dimension::T6
    "Dimension of shared inequality constraint"
    shared_inequality_dimension::T6

    "Corresponding ParametricMCP."
    parametric_mcp::ParametricMCP
end

function ParametricGame(;
    objectives,
    equality_constraints,
    inequality_constraints,
    shared_equality_constraint,
    shared_inequality_constraint,
    parameter_dimension = 1,
    primal_dimensions,
    equality_dimensions,
    inequality_dimensions,
    shared_equality_dimension,
    shared_inequality_dimension,
)
    @assert !isnothing(equality_constraints)
    @assert !isnothing(inequality_constraints)

    N = length(objectives)
    @assert N ==
            length(equality_constraints) ==
            length(inequality_constraints) ==
            length(primal_dimensions) ==
            length(equality_dimensions) ==
            length(inequality_dimensions)

    total_dimension =
        sum(primal_dimensions) +
        sum(equality_dimensions) +
        sum(inequality_dimensions) +
        shared_equality_dimension +
        shared_inequality_dimension

    # Define symbolic variables for this MCP.
    # TODO!
    z̃ = Symbolics.scalarize(only(@variables(z̃[1:total_dimension])))
    z = BlockArray(z̃, [sum(primal_dimensions), sum(equality_dimensions), sum(inequality_dimensions), 
                       shared_equality_dimension, shared_inequality_dimension])
    
    x = BlockArray(z[Block(1)], primal_dimensions)
    λ = BlockArray(z[Block(2)], equality_dimensions)
    μ = BlockArray(z[Block(3)], inequality_dimensions)
    λ̃ = z[Block(4)]
    μ̃ = z[Block(5)]

    # Define a symbolic variable for the parameters.
    # TODO!
    @variables θ̃[1:(parameter_dimension)]
    θ = Symbolics.scalarize(θ̃)

    # Build symbolic expressions for objectives and constraints for all players
    # (and shared constraints).
    # TODO!
    fs = [obj for obj in objectives]
    gs = [eq for eq in equality_constraints]
    hs = [ineq for ineq in inequality_constraints]

    # handle shared constraints like in parametric_optimization_problem.jl
    g̃ = shared_equality_constraint(x, θ)
    h̃ = shared_inequality_constraint(x, θ)

    # Build Lagrangians for all players.
    # TODO!
    grads = []
    for i in 1:N
        equality = λ[Block(i)]' * gs[i](x, θ)
        inequality = μ[Block(i)]' * hs[i](x, θ)
        L = fs[i](x, θ) - equality - inequality - λ̃' * g̃ - μ̃' * h̃
        grad = Symbolics.gradient(L, x[Block(i)])
        push!(grads, grad)
    end

    # Build F = [∇ₓLs, gs, hs, g̃, h̃]'.
    # F = # TODO!
    @variables Fvec[1:total_dimension]
    Fvec = Symbolics.scalarize(Fvec)

    # add components to F
    F = []
    [push!(F, val) for grad in grads for val in grad]
    [push!(F, val) for g in gs for val in g(x, θ)]
    [push!(F, val) for h in hs for val in h(x, θ)]
    [push!(F, val) for val in g̃]
    [push!(F, val) for val in h̃]

    # fill in Fvec from F
    [Fvec[i] = F[i] for i in 1:length(F)]

    # Set lower and upper bounds for z.
    # z̲ = # TODO!
    # z̅ = # TODO!
    z̲ = Vector{Float64}(undef, 0)
    z̅ = Vector{Float64}(undef, 0)

    # all upper bounds are inf
    [push!(z̅, Inf) for _ in 1:total_dimension]

    # lower bounds are zero for inequality constraints, -inf for everything else
    [push!(z̲, -Inf) for _ in 1:sum(primal_dimensions)]
    [push!(z̲, -Inf) for _ in 1:sum(equality_dimensions)]
    [push!(z̲, 0) for _ in 1:sum(inequality_dimensions)]
    [push!(z̲, -Inf) for _ in 1:shared_equality_dimension]
    [push!(z̲, 0) for _ in 1:shared_inequality_dimension]

    # Build parametric MCP.
    # parametric_mcp = ParametricMCP(F, z̲, z̅, parameter_dimension)
    parametric_mcp = ParametricMCP(Fvec, z̃, θ, z̲, z̅; compute_sensitivities = false)

    ParametricGame(
        objectives,
        equality_constraints,
        inequality_constraints,
        shared_equality_constraint,
        shared_inequality_constraint,
        parameter_dimension,
        primal_dimensions,
        equality_dimensions,
        inequality_dimensions,
        shared_equality_dimension,
        shared_inequality_dimension,
        parametric_mcp,
    )
end

function total_dim(problem::ParametricGame)
    sum(problem.primal_dimensions) +
    sum(problem.equality_dimensions) +
    sum(problem.inequality_dimensions) +
    problem.shared_equality_dimension +
    problem.shared_inequality_dimension
end

"Solve a constrained parametric game."
function solveGame(
    problem::ParametricGame,
    parameter_value = zeros(problem.parameter_dimension);
    initial_guess = nothing,
    verbose = false,
)
    z0 = if !isnothing(initial_guess)
        initial_guess
    else
        zeros(total_dim(problem))
    end

    z, status, info = ParametricMCPs.solve(
        problem.parametric_mcp,
        parameter_value;
        initial_guess = z0,
        verbose,
        cumulative_iteration_limit = 100000,
        proximal_perturbation = 1e-2,
        use_basics = true,
        use_start = true,
    )

    primals = blocks(BlockArray(z[1:sum(problem.primal_dimensions)], problem.primal_dimensions))
    (; primals, variables = z, status, info)
end
""" Write a standard-form constrained static game as a MCP
i.e., if player i's problem is of the form
             min_{xᵢ} fᵢ(x, θ)
             s.t.     gᵢ(x, θ) = 0
                      hᵢ(x, θ) ≥ 0

and players share constraints
                      g̃(x, θ) = 0
                      h̃(x, θ) ≥ 0

Here, xᵢ is Pi's decision variable, x is the concatenation of all players'
variables, and θ is a vector of parameters. Express it in the following form:
             find   z
             s.t.   F(z, θ) ⟂ z̲ ≤ z ≤ z̅

where we interpret z = (x, λ, μ, λ̃, μ̃), with λ and μ the Lagrange multipliers
for the constraints g and h (with tildes as appropriate), respectively. The
expression F(z) ⟂ z̲ ≤ z ≤ z̅ should be read as the following three statements:
            - if z = z̲, then F(z, θ) ≥ 0
            - if z̲ < z < z̅, then F(z, θ) = 0
            - if z = z̅, then F(z, θ) ≤ 0

For more details, please consult the documentation for the package
`Complementarity.jl`, which may be found here:
https://github.com/chkwon/Complementarity.jl/tree/master
"""

"Generic description of a constrained parametric game problem."
struct ParametricGame{T1,T2,T3,T4,T5,T6,T7}
    "Objective functions for all players"
    objectives::T1
    "Equality constraints for all players"
    equality_constraints::T2
    "Inequality constraints for all players"
    inequality_constraints::T3
    "Shared equality constraint"
    shared_equality_constraint::T4
    "Shared inequality constraint"
    shared_inequality_constraint::T5

    "Dimension of parameter vector"
    parameter_dimension::T6
    "Dimension of primal variables for all players"
    primal_dimensions::T7
    "Dimension of equality constraints for all players"
    equality_dimensions::T7
    "Dimension of inequality constraints for all players"
    inequality_dimensions::T7
    "Dimension of shared equality constraint"
    shared_equality_dimension::T6
    "Dimension of shared inequality constraint"
    shared_inequality_dimension::T6

    "Corresponding ParametricMCP."
    parametric_mcp::ParametricMCP
end

function ParametricGame(;
    objectives,
    equality_constraints,
    inequality_constraints,
    shared_equality_constraint,
    shared_inequality_constraint,
    parameter_dimension = 1,
    primal_dimensions,
    equality_dimensions,
    inequality_dimensions,
    shared_equality_dimension,
    shared_inequality_dimension,
)
    @assert !isnothing(equality_constraints)
    @assert !isnothing(inequality_constraints)

    N = length(objectives)
    @assert N ==
            length(equality_constraints) ==
            length(inequality_constraints) ==
            length(primal_dimensions) ==
            length(equality_dimensions) ==
            length(inequality_dimensions)

    total_dimension =
        sum(primal_dimensions) +
        sum(equality_dimensions) +
        sum(inequality_dimensions) +
        shared_equality_dimension +
        shared_inequality_dimension

    # Define symbolic variables for this MCP.
    # TODO!
    z̃ = Symbolics.scalarize(only(@variables(z̃[1:total_dimension])))
    z = BlockArray(z̃, [sum(primal_dimensions), sum(equality_dimensions), sum(inequality_dimensions), 
                       shared_equality_dimension, shared_inequality_dimension])
    
    x = BlockArray(z[Block(1)], primal_dimensions)
    λ = BlockArray(z[Block(2)], equality_dimensions)
    μ = BlockArray(z[Block(3)], inequality_dimensions)
    λ̃ = z[Block(4)]
    μ̃ = z[Block(5)]

    # Define a symbolic variable for the parameters.
    # TODO!
    @variables θ̃[1:(parameter_dimension)]
    θ = Symbolics.scalarize(θ̃)

    # Build symbolic expressions for objectives and constraints for all players
    # (and shared constraints).
    # TODO!
    fs = [obj for obj in objectives]
    gs = [eq for eq in equality_constraints]
    hs = [ineq for ineq in inequality_constraints]

    # handle shared constraints like in parametric_optimization_problem.jl
    g̃ = shared_equality_constraint(x, θ)
    h̃ = shared_inequality_constraint(x, θ)

    # Build Lagrangians for all players.
    # TODO!
    grads = []
    for i in 1:N
        equality = λ[Block(i)]' * gs[i](x, θ)
        inequality = μ[Block(i)]' * hs[i](x, θ)
        L = fs[i](x, θ) - equality - inequality - λ̃' * g̃ - μ̃' * h̃
        grad = Symbolics.gradient(L, x[Block(i)])
        push!(grads, grad)
    end

    # Build F = [∇ₓLs, gs, hs, g̃, h̃]'.
    # F = # TODO!
    @variables Fvec[1:total_dimension]
    Fvec = Symbolics.scalarize(Fvec)

    # add components to F
    F = []
    [push!(F, val) for grad in grads for val in grad]
    [push!(F, val) for g in gs for val in g(x, θ)]
    [push!(F, val) for h in hs for val in h(x, θ)]
    [push!(F, val) for val in g̃]
    [push!(F, val) for val in h̃]

    # fill in Fvec from F
    [Fvec[i] = F[i] for i in 1:length(F)]

    # Set lower and upper bounds for z.
    # z̲ = # TODO!
    # z̅ = # TODO!
    z̲ = Vector{Float64}(undef, 0)
    z̅ = Vector{Float64}(undef, 0)

    # all upper bounds are inf
    [push!(z̅, Inf) for _ in 1:total_dimension]

    # lower bounds are zero for inequality constraints, -inf for everything else
    [push!(z̲, -Inf) for _ in 1:sum(primal_dimensions)]
    [push!(z̲, -Inf) for _ in 1:sum(equality_dimensions)]
    [push!(z̲, 0) for _ in 1:sum(inequality_dimensions)]
    [push!(z̲, -Inf) for _ in 1:shared_equality_dimension]
    [push!(z̲, 0) for _ in 1:shared_inequality_dimension]

    # Build parametric MCP.
    # parametric_mcp = ParametricMCP(F, z̲, z̅, parameter_dimension)
    parametric_mcp = ParametricMCP(Fvec, z̃, θ, z̲, z̅; compute_sensitivities = false)

    ParametricGame(
        objectives,
        equality_constraints,
        inequality_constraints,
        shared_equality_constraint,
        shared_inequality_constraint,
        parameter_dimension,
        primal_dimensions,
        equality_dimensions,
        inequality_dimensions,
        shared_equality_dimension,
        shared_inequality_dimension,
        parametric_mcp,
    )
end

function total_dim(problem::ParametricGame)
    sum(problem.primal_dimensions) +
    sum(problem.equality_dimensions) +
    sum(problem.inequality_dimensions) +
    problem.shared_equality_dimension +
    problem.shared_inequality_dimension
end

"Solve a constrained parametric game."
function solve(
    problem::ParametricGame,
    parameter_value = zeros(problem.parameter_dimension);
    initial_guess = nothing,
    verbose = false,
)
    z0 = if !isnothing(initial_guess)
        initial_guess
    else
        zeros(total_dim(problem))
    end

    z, status, info = ParametricMCPs.solve(
        problem.parametric_mcp,
        parameter_value;
        initial_guess = z0,
        verbose,
        cumulative_iteration_limit = 100000,
        proximal_perturbation = 1e-2,
        use_basics = true,
        use_start = true,
    )

    primals = blocks(BlockArray(z[1:sum(problem.primal_dimensions)], problem.primal_dimensions))
    (; primals, variables = z, status, info)
end
