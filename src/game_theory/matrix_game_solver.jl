abstract type FiniteGameSolver end

## ====================================================================

abstract type AbstractTrajectoryGenerator end
export AbstractTrajectoryGenerator

## ====================================================================

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

## ====================================================================

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
    (; x=sol1.x, y=sol2.x)
end

export solve_mixed_nash

## ====================================================================

function solve_mixed_security_strategy(player_cost_matrix)

    # TODO: transform the game to ensure that the cost matrix is entrywise positive
    r = size(player_cost_matrix, 1)
    p = size(player_cost_matrix, 2)

    min_value = 0
    for i in 1:r
        for j in 1:p
            if player_cost_matrix[i, j] <= min_value
                min_value = player_cost_matrix[i, j]
            end
        end
    end
    if min_value <= 0
        c = -min_value + 1
        M = player_cost_matrix + (-min_value + 1) * ones(r, p)
    end

    # TODO: solve the LP associated to a zero sum game
    ans = solve_simplex_lp(M)
    x_tilde = ans.x
    V_tilde = ans.V

    # TODO: transform the solution into the probability simplex
    x_star = x_tilde * (1 / V_tilde)
    V_star = (1 / V_tilde) - c

    # @infiltrate

    # TODO: return a named tuple of (; x, V) where x is the strategy and V is the value
    (; x=x_star, v=V_star)
end

export solve_mixed_security_strategy

## ====================================================================

function solve_simplex_lp(A)

    # set-up the optimization problem 
    model = JuMP.Model()
    JuMP.set_optimizer(model, OSQP.Optimizer)
    JuMP.set_silent(model)

    # get dimensions 
    r, p = size(A)

    # TODO: add constraints and objective
    @variable(model, z[1:r])
    @objective(model, Max, ones(r)' * z)
    @constraint(model, c1, ones(p) >= A' * z)

    for i in 1:r
        @constraint(model, z[i] >= 1e-4)
    end

    JuMP.optimize!(model)

    (JuMP.termination_status(model) == JuMP.MOI.OPTIMAL) ||
    # error("OSQP did not find an optimal solution to this matrix game.")
        println("termination status = ", JuMP.termination_status(model))
    (; x=JuMP.value.(z), V=JuMP.objective_value(model))
end

export solve_simplex_lp

## ====================================================================

# compute X, U, and t for a player (optimization) 
function compute_players_XU(params, game, players)

    # get vertices 
    rv_ref_polygon = game.rv_ref_E[end][end, :]
    vertices = polygon_vertices(rv_ref_polygon, params)

    for ii in eachindex(players)

        p = players[ii]
        rv_0_hist = p.rv_0_hist

        # init state and end velocity (probably doesn't matter) 
        rv_0 = rv_0_hist[1, :]
        v_f = rv_0_hist[end, 4:6]

        # get params 
        tof = params.tof
        N = params.N
        mu = params.mu

        for jj in eachindex(vertices)

            rv_f = [vertices[jj]; v_f]
            Δv_sol = min_Δv_dist(rv_0, rv_f, tof, N, mu)
            t, rv_hist = prop_kepler_tof_Nseg(rv_0, Δv_sol, N, tof / N, mu)

            # save hist 
            push!(p.X, rv_hist)
            push!(p.U, Δv_sol)
            push!(p.t, t)

        end

    end

    return players
end

export compute_players_XU


## ====================================================================

# game cost 
function stage_cost(x1, x2, u1, u2)

    # sqrt( norm(x1[1:3] - x2[1:3]) + 0.1 ) + 0.1 * (norm(u1) - norm(u2))
    dist = norm(x1[1:3] - x2[1:3])
    cost = 1.0 * sqrt(dist + 0.1) + 0.1 * (norm(u1 - u2))

    # if dist < capture_threshold
    #     cost -= 50.0
    # end

    return cost
end

export stage_cost


## ==================================================================== 
# zero-sum game 

function compute_cost_matrices(players)

    # loop through time corresponding with control inputs 
    U_idx = eachindex(players[1].U[1][:, 1])

    # start with vertex 1 for players 1 and 2 
    i_vert = 1
    j_vert = 1

    player1_cost_matrix = zeros(6, 6)
    player2_cost_matrix = zeros(6, 6)

    n_vertices = length(players[begin].t)

    for i_vert in 1:n_vertices
        for j_vert in 1:n_vertices

            # loop through time 
            player1_cost_tt = []
            player2_cost_tt = []

            for ii in U_idx

                x1 = players[1].X[i_vert][ii, :]
                u1 = players[1].U[i_vert][ii, :]
                x2 = players[2].X[j_vert][ii, :]
                u2 = players[2].U[j_vert][ii, :]

                # compute costs for player 1 and 2  
                cost1 = stage_cost(x1, x2, u1, u2)
                cost2 = -stage_cost(x1, x2, u1, u2)
                push!(player1_cost_tt, cost1)
                push!(player2_cost_tt, cost2)

            end

            # @infiltrate

            player1_cost = mean(player1_cost_tt)
            player2_cost = mean(player2_cost_tt)
            # player1_cost = player1_cost_tt[end]
            # player2_cost = player2_cost_tt[end]

            # save cost in matrix 
            player1_cost_matrix[i_vert, j_vert] = player1_cost
            player2_cost_matrix[i_vert, j_vert] = player2_cost

        end
    end

    players[1].cost = player1_cost_matrix
    players[2].cost = player2_cost_matrix

    return players
end

export compute_cost_matrices


## ====================================================================

function stage_cost_games_fn(games_vec)

    stage_cost_games = []
    # stage_cost_games_mean = [ ]
    for ii in eachindex(games_vec)

        game = games_vec[ii]
        params = game.params[1]

        p1_U_hist, p2_U_hist = p1_p2_u_hist(game)
        tt_hist, p1_rv_hist, p2_rv_hist, rv_ref_hist = p_rv_ref_hist(game, params)

        x1 = p1_rv_hist
        x2 = p2_rv_hist
        u1 = p1_U_hist
        u2 = p2_U_hist

        # compute stage cost at each time step for each trajectory 
        stage_cost_game = [stage_cost(x1[ii, :], x2[ii, :], u1[ii, :], u2[ii, :]) for ii in 1:size(u1, 1)]
        # stage_cost_game_mean = mean( stage_cost_game ) 

        push!(stage_cost_games, stage_cost_game)
        # push!( stage_cost_games_mean, stage_cost_game_mean ) 

    end

    stage_cost_games = mapreduce(permutedims, vcat, stage_cost_games)
    stage_cost_games_mean = mean(stage_cost_games, dims=1)[:]
    stage_cost_games_std = std(stage_cost_games, dims=1)[:]

    return stage_cost_games, stage_cost_games_mean, stage_cost_games_std
end

export stage_cost_games_fn


## ==================================================================== 

# using Infiltrator

function compute_mixing_weights!(players)

    # mixing weights - ZERO SUM GAME!!! 
    mixing_weights = let
        sol = solve_mixed_nash(players[1].cost)
        (; sol.x, sol.y)
    end
    players[1].weights = mixing_weights[1]
    players[2].weights = mixing_weights[2]

    return mixing_weights
end

## ====================================================================

function choose_strategies!(players, mixing_weights, rng, params, game)

    # now determine strategy 
    p1_strategy = params.strategy
    p2_strategy = params.p2_strategy

    # initialize chosen vector 
    chosen = [0, 0]

    # ---------------------------------- 
    # determine strategy for player 1 
    # ---------------------------------- 

    if p1_strategy == "mixed"
        chosen[1] = sample(rng, ProbabilityWeights(mixing_weights[1]))

    elseif p1_strategy == "greedy"
        chosen[1] = argmax(mixing_weights[1])

    elseif p1_strategy == "random"
        len = length(mixing_weights[1])
        chosen[1] = rand(rng, 1:len)

    elseif p1_strategy == "FP_greedy"
        # P1's belief about P2's moves → best response 
        q_prob = players[1].fp_belief ./ sum(players[1].fp_belief)
        expected_costs = players[1].cost * q_prob
        chosen[1] = argmax(expected_costs)

    elseif p1_strategy == "FP_mixed"
        # convert to weights and ensure positive 
        q_prob = players[1].fp_belief ./ sum(players[1].fp_belief)
        expected_costs = players[1].cost * q_prob
        w = expected_costs .- minimum(expected_costs) .+ 1e-6
        chosen[1] = sample(rng, ProbabilityWeights(w))

    elseif p1_strategy == "Meta_greedy"
        predicted_v_probs = predict_opponent_vertices(players[1], players[2])
        expected_costs = players[1].cost * predicted_v_probs
        chosen[1] = argmax(expected_costs)

    elseif p1_strategy == "Meta_mixed"
        predicted_v_probs = predict_opponent_vertices(players[1], players[2])
        meta_expected_costs = players[1].cost * predicted_v_probs
        w = meta_expected_costs .- minimum(meta_expected_costs) .+ 1e-6
        chosen[1] = sample(rng, ProbabilityWeights(w))

    elseif p1_strategy isa Int
        chosen[1] = p1_strategy

    else
        error("Invalid p1 strategy: $p1_strategy")
    end

    # ---------------------------------- 
    # determine strategy for player 2 
    # ---------------------------------- 

    if p2_strategy == "mixed"
        chosen[2] = sample(rng, ProbabilityWeights(mixing_weights[2]))

    elseif p2_strategy == "greedy"
        chosen[2] = argmax(mixing_weights[2])

    elseif p2_strategy == "random"
        len = length(mixing_weights[2])
        chosen[2] = rand(rng, 1:len)

    elseif p2_strategy == "FP_greedy"
        q_prob = players[2].fp_belief / sum(players[2].fp_belief)
        expected_costs = players[2].cost' * q_prob
        chosen[2] = argmax(expected_costs)

    elseif p2_strategy == "FP_mixed"
        # convert to weights and ensure positive 
        q_prob = players[2].fp_belief / sum(players[2].fp_belief)
        expected_costs = players[2].cost' * q_prob
        w = expected_costs .- minimum(expected_costs) .+ 1e-6
        chosen[2] = sample(rng, ProbabilityWeights(w))

    elseif p2_strategy == "Meta_greedy"
        predicted_v_probs = predict_opponent_vertices(players[2], players[1])
        expected_costs = players[2].cost * predicted_v_probs
        # player 2 minimizes so we take argmin 
        chosen[2] = argmin(expected_costs)

    elseif p2_strategy == "Meta_mixed"
        predicted_v_probs = predict_opponent_vertices(players[2], players[1])
        meta_expected_costs = players[2].cost * predicted_v_probs
        # for player 2 we negate or max over negative Expected Costs
        w = maximum(meta_expected_costs) .- meta_expected_costs .+ 1e-6
        chosen[2] = sample(rng, ProbabilityWeights(w))

    elseif p2_strategy isa Int
        chosen[2] = p2_strategy

    else
        error("Invalid p2 strategy: $p2_strategy")
    end

    # ---------------------------------- 
    # save chosen vertex 
    # ---------------------------------- 

    players[1].chosen = chosen[1]
    players[2].chosen = chosen[2]

    return chosen
end

## ====================================================================

function update_chosen_trajectories!(players, params)

    for i in 1:2
        idx = players[i].chosen
        players[i].t_chosen = players[i].t[idx][1:params.k_tt_replan+1]
        players[i].rv_chosen = players[i].X[idx][1:params.k_tt_replan+1, :]
        players[i].U_chosen = players[i].U[idx][1:params.k_tt_replan, :]
    end

    return players
end

## ====================================================================

"""
Update fictitious play belief based on the opponent's chosen vertex.
Increments the count for the vertex the opponent actually chose.
"""
function update_beliefs!(game, players)

    # P1 observes P2's choice → update P1's belief about P2
    players[1].fp_belief[players[2].chosen] += 1
    # P2 observes P1's choice → update P2's belief about P1
    players[2].fp_belief[players[1].chosen] += 1

    # Bayesian update of meta-strategy belief
    update_strategy_belief!(players[1], players[2])
    update_strategy_belief!(players[2], players[1])

end

## ====================================================================

function update_strategy_belief!(p_self, p_opponent)

    # 1. CALCULATING THE LIKELIHOODS: P(v_t | s)
    # How likely was the opponent's move if they were using strategy 's'?
    likelihoods = zeros(length(p_self.tracked_strategies))

    for (i, s) in enumerate(p_self.tracked_strategies)
        if s == "mixed"
            # P(v_t | s) = w[v_t]
            likelihoods[i] = p_opponent.weights[p_opponent.chosen]

        elseif s == "greedy"
            # P(v_t | s) = 1 if max weight, else 0
            likelihoods[i] = (p_opponent.chosen == argmax(p_opponent.weights)) ? 1.0 : 0.0

        elseif s == "random"
            # P(v_t | s) = 1 / N_v
            likelihoods[i] = 1.0 / length(p_opponent.weights)

        elseif s isa Int
            # P(v_t | s) = 1 if v_t == k, else 0
            likelihoods[i] = (p_opponent.chosen == s) ? 1.0 : 0.0

        else
            likelihoods[i] = 1.0 / length(p_opponent.weights) # fallback
        end
    end

    # 2. ADDING EPSILON
    # Math note: Add 10^-4 so beliefs never hit absolute zero.
    likelihoods = max.(likelihoods, 1e-4)

    # 3. APPLYING BAYES' THEOREM: P(S|v) = P(v|S) * P(S) / sum(...)

    # Numerator: P(v_t | s) * P_t(s)
    # Multiplies our current belief by the likelihood we just calculated.
    p_self.strategy_belief .*= likelihoods

    # Denominator (Normalization): Divide by the sum of all probabilities
    # This ensures all our new belief probabilities still add up to 1.0 (100%).
    p_self.strategy_belief ./= sum(p_self.strategy_belief)

end

## ====================================================================

function predict_opponent_vertices(p_self, p_opponent)

    n_v = length(p_opponent.weights)
    predicted_v_probs = zeros(n_v) # This will be our final 'p' vector

    for (i, s) in enumerate(p_self.tracked_strategies)

        # P_{t+1}(s) 
        # Grabbing our updated confidence in strategy 's' from Step 1.
        prob_s = p_self.strategy_belief[i]

        v_probs = zeros(n_v)

        # P(v_{t+1} = i | s)
        # Figuring out the opponent's vertex probabilities IF they are using 's'
        if s == "mixed"
            v_probs .= p_opponent.weights
        elseif s == "greedy"
            v_probs[argmax(p_opponent.weights)] = 1.0
        elseif s == "random"
            v_probs .= 1.0 / n_v
        elseif s isa Int
            if s >= 1 && s <= n_v
                v_probs[s] = 1.0
            end
        else
            v_probs .= 1.0 / n_v # fallback
        end

        # THE CORE EQUATION: Summation part!
        # Math: p_i += P(v_{t+1} = i | s) * P_{t+1}(s)
        # We multiply their hypothetical move (v_probs) by our belief in that 
        # hypothesis (prob_s), and add it to our running total.
        predicted_v_probs .+= prob_s .* v_probs
    end

    # Returns the final 'p' vector (normalized just to be safe)
    return predicted_v_probs / sum(predicted_v_probs)
end

export compute_mixing_weights!, choose_strategies!, update_chosen_trajectories!, update_beliefs!, predict_opponent_vertices


## ====================================================================

# compute X, U, and t for a player 
function compute_states_nash(params, game, players, rng)

    # compute X, U, and t for both players (optimization)   
    players = compute_players_XU(params, game, players)

    # compute cost matrices 
    players = compute_cost_matrices(players)

    # solve mixed nash and choose strategy 
    mixing_weights = compute_mixing_weights!(players)

    # choose strategies 
    choose_strategies!(players, mixing_weights, rng, params, game)

    # update beliefs 
    update_beliefs!(game, players)

    # update chosen trajectories  
    players = update_chosen_trajectories!(players, params)

    return players
end

export compute_states_nash

