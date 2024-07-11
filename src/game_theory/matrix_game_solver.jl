abstract type FiniteGameSolver end

## ============================================ ##

"""
    solve_mixed_nash(solver, A)
The entry-point to a game solver.
Inputs:
- solver: low-level solver to be used (e.g. MatrixGameSolver)
- A: the cost matrix for P1.
Returns a named tuple of
- x: strategy for P1
- y: strategy for P2
- V: the Nash value (for P1)
"""
function solve_mixed_nash end

"""
Returns the game cost associated to strategies x and y for P1 and P2, respectively.
A is the cost matrix for P1.
"""
function game_cost(x, y, A)
    x' * A * y
end

"""
A zero-sum game solver that casts the game as linear program.
"""
struct MatrixGameSolver <: FiniteGameSolver end

# function solve_mixed_nash(A, B)
#     solve_mixed_nash( A )
# end 

function solve_mixed_nash(A)
    sol1 = solve_mixed_security_strategy(A)
    sol2 = solve_mixed_security_strategy(-A')
    (; x = sol1.x, y = sol2.x)
end

export solve_mixed_nash 

## ============================================ ##

function solve_mixed_security_strategy(player_cost_matrix)

    # TODO: transform the game to ensure that the cost matrix is entrywise positive
    r = size(player_cost_matrix, 1)
    p = size(player_cost_matrix, 2)
    min_value = 0 
    for i in 1 : r
        for j in 1 : p
            if player_cost_matrix[i,j] <= min_value
                min_value = player_cost_matrix[i,j]
            end
        end
    end
    if min_value <= 0
        c = -min_value+1
        M = player_cost_matrix + (-min_value+1) * ones(r,p) 
    end

    # TODO: solve the LP associated to a zero sum game
    ans = solve_simplex_lp(M)
    x_tilde = ans.x
    V_tilde = ans.V
    
    # TODO: transform the solution into the probability simplex
    x_star = x_tilde * (1 / V_tilde)
    V_star = (1 / V_tilde) - c

    # TODO: return a named tuple of (; x, V) where x is the strategy and V is the value
    (; x = x_star, v = V_star)
end 

export solve_mixed_security_strategy 

## ============================================ ##

function solve_simplex_lp(A)

    # set-up the optimization problem 
    model = JuMP.Model()
    JuMP.set_optimizer(model, OSQP.Optimizer)
    JuMP.set_silent(model)

    # get dimensions 
    r, p = size(A)

    # TODO: add constraints and objective
    @variable( model, z[1 : r] )
    @objective( model, Max, ones(r)' * z ) 
    @constraint( model, c1, ones(p) >= A'*z )

    for i in 1 : r 
        @constraint( model,z[i] >= 1e-4 )
    end

    JuMP.optimize!(model)
    (JuMP.termination_status(model) == JuMP.MOI.OPTIMAL) ||
        error("OSQP did not find an optimal solution to this matrix game.")
    (; x = JuMP.value.(z), V = JuMP.objective_value(model))
end

export solve_simplex_lp 
