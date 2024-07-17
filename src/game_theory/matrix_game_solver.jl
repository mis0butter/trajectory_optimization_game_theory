abstract type FiniteGameSolver end

## ============================================ ##

abstract type AbstractTrajectoryGenerator end
export AbstractTrajectoryGenerator 

## ============================================ ##

struct FiniteTrajectoryGameSolver{TG,TT<:AbstractTrajectoryGenerator,TH,TR,TF<:FiniteGameSolver}
    "Underlying trajectory game to be solved."
    game::TG
    "A trajectory generator to be used by all players."
    trajectory_generator::TT
    "The number of time steps to plan into the future."
    planning_horizon::TH
    "A random number generator to generate non-deterministic strategies."
    rng::TR
    "The solver for the high-level finite game."
    finite_game_solver::TF
end

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

## ============================================ ##

# compute X, U, and t for a player 
function player_XU( params, vertices, players ) 

    for ii in eachindex(players) 

        p = players[ii] 
        rv_player = p.rv_0 

        # init state and end velocity (probably doesn't matter) 
        rv_0 = rv_player[1,:] 
        v_f  = rv_player[end,4:6] 
    
        # get params 
        tof = params.tof ; N = params.N ; mu  = params.mu 
    
        for jj in eachindex(vertices)  
    
            rv_f   = [ vertices[jj] ; v_f ]  
            Δv_sol = min_Δv_dist( rv_0, rv_f, tof, N, mu ) 
            t, rv_hist = prop_kepler_tof_Nseg( rv_0, Δv_sol, N, tof / N, mu ) 
    
            # save hist 
            push!(p.X, rv_hist) 
            push!(p.U, Δv_sol) 
            push!(p.t, t) 
    
        end 

    end 

    return players 
end 

export player_XU 


## ============================================ ## 
# zero-sum game 

function player_cost_matrices( players ) 

    # game cost 
    function stage_cost(x1, x2, u1, u2)
        sqrt(norm(x1[1:3] - x2[1:3]) + 0.1) + 0.1 * (norm(u1) - norm(u2))
    end

    # loop through time corresponding with control inputs 
    U_idx = eachindex( players[1].U[1][:,1] )

    # start with vertex 1 for players 1 and 2 
    i_vert = 1 
    j_vert = 1 

    player1_cost_matrix = zeros(6, 6) 
    player2_cost_matrix = zeros(6, 6) 

    n_vertices = length( players[begin].t ) 

    for i_vert in 1 : n_vertices 
        for j_vert in 1 : n_vertices 

            # loop through time 
            player1_cost_tt = [] 
            player2_cost_tt = [] 
            for ii in U_idx 

                x1 = players[1].X[i_vert][ii,:] 
                u1 = players[1].U[i_vert][ii,:] 
                x2 = players[2].X[j_vert][ii,:] 
                u2 = players[2].U[j_vert][ii,:] 

                # compute costs for player 1 and 2  
                cost1 = stage_cost( x1, x2, u1, u2 )
                cost2 = - stage_cost( x1, x2, u1, u2 )
                push!( player1_cost_tt, cost1 ) 
                push!( player2_cost_tt, cost2 ) 
                
            end 
            player1_cost = mean( player1_cost_tt )
            player2_cost = mean( player2_cost_tt )  

            # save cost in matrix 
            player1_cost_matrix[i_vert, j_vert] = player1_cost 
            player2_cost_matrix[i_vert, j_vert] = player2_cost 

        end 
    end 

    players[1].cost = player1_cost_matrix   
    players[2].cost = player2_cost_matrix   

    return players 
end 

export player_cost_matrices 


## ============================================ ##

function choose_weights( players, rng ) 

    # mixing weights - ZERO SUM GAME!!! 
    mixing_weights = let
        sol = solve_mixed_nash( players[1].cost ) 
        (; sol.x, sol.y) 
    end 
    players[1].weights = mixing_weights[1] 
    players[2].weights = mixing_weights[2] 

    # sample from mixed nash 
    chosen = [sample(rng, ProbabilityWeights(weights)) for weights in mixing_weights] 
    players[1].chosen = chosen[1] 
    players[2].chosen = chosen[2] 

    return players 
end 

export choose_weights 



