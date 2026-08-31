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

export game_cost

"""
A zero-sum game solver that casts the game as linear program.
"""
struct MatrixGameSolver <: FiniteGameSolver end

# function solve_mixed_nash(A, B)
#     solve_mixed_nash( A )
# end 

"Project a raw LP iterate onto the probability simplex (clamp, then normalize)."
function normalize_simplex(w)
    v = max.(vec(w), 0.0)
    s = sum(v)
    s > 1e-12 ? v ./ s : fill(inv(length(v)), length(v))
end

export normalize_simplex

"""
    nash_certificate(A, x, y, V)

Solver-independent check that `(x, y, V)` really is a saddle point of the
zero-sum game with P1 payoff `A`, where P1 maximizes and P2 minimizes.

Deliberately consults no solver status flag. OSQP is a first-order ADMM method
being used on a pure LP; asking it whether it converged is weaker evidence than
checking the optimality conditions directly. Reported as a scale-free residual so
it can be aggregated across game steps with wildly different cost magnitudes.
"""
function nash_certificate(A, x, y, V)
    scale = max(1.0, maximum(abs, A))
    v_bilinear = game_cost(x, y, A)
    violation = max(
        abs(sum(x) - 1),                # x on the simplex
        abs(sum(y) - 1),                # y on the simplex
        max(0.0, -minimum(x)),          # x nonnegative
        max(0.0, -minimum(y)),          # y nonnegative
        max(0.0, maximum(A * y) - V),   # y holds P1 to at most V
        max(0.0, V - minimum(A' * x)),  # x guarantees P1 at least V
        abs(v_bilinear - V),            # value agrees with the bilinear form
    ) / scale
    (; violation, v_bilinear)
end

export nash_certificate

"""
    solve_mixed_nash(A)

Zero-sum matrix game on P1's payoff matrix `A`: `A[i,j]` is P1's payoff when P1
plays row `i` and P2 plays column `j`.

**Orientation (defect D1).** `stage_cost` increases with the separation between
the two spacecraft, and `compute_cost_matrices` sets `players[1].cost = +stage_cost`,
`players[2].cost = -stage_cost`. Player 1 is the evader, who wants separation, so
**P1 MAXIMIZES `A` and P2 MINIMIZES `A`** — which is exactly what every
fictitious-play branch in `choose_strategies!` already assumes (`argmax` on P1's
expected costs).

`solve_mixed_security_strategy(M)` solves `min_x max_j (x'M)_j`: it returns the
security strategy of the row player who *minimizes* `M`. Hence

    P1 maximizes A  <=>  P1 minimizes -A   ->  sms(-A),  maximin = -v
    P2 minimizes A over its columns        ->  sms(A'),  minimax = +v

The previous implementation called `sms(A)` and `sms(-A')` — both of these with
the orientation inverted, handing the evader a distance-*minimizing* mixed
strategy and the pursuer a distance-*maximizing* one. Every `mixed` and `greedy`
number in the CDC submission is therefore the equilibrium of the reversed game.
Reviewer 7 (#2, #3), Reviewer 8 (#1) and the AE all flagged this.

By von Neumann's minimax theorem maximin == minimax; `gap` is the numerical
evidence, and `cert_violation` independently certifies the saddle point.
"""
function solve_mixed_nash(A)

    sol1 = solve_mixed_security_strategy(-A)   # P1: the MAXIMIZER of A
    sol2 = solve_mixed_security_strategy(A')   # P2: the MINIMIZER of A

    x, y      = normalize_simplex(sol1.x), normalize_simplex(sol2.x)
    V_lo, V_hi = -sol1.v, sol2.v               # maximin, minimax
    V          = 0.5 * (V_lo + V_hi)

    cert = nash_certificate(A, x, y, V)

    # Failure policy: never throw. This runs inside `Threads.@threads` in
    # run_MC_games_parallel, where an error would take down a whole sweep rather
    # than a single 6x6 solve. Instead, fall back to the pure maximin/minimax
    # strategies -- always well defined, deterministic -- and mark the record so
    # the failure rate can be reported rather than silently absorbed.
    fallback = false
    if cert.violation > 1e-8
        # Pure maximin / minimax, computed directly from A so they are exact
        # regardless of what the solver did. These bracket the true value:
        # V_lo_pure <= V <= V_hi_pure, so the midpoint is a principled estimate,
        # unlike the single matrix entry game_cost(x,y,A) would give.
        row_min, col_max = vec(minimum(A, dims=2)), vec(maximum(A, dims=1))
        x = zeros(size(A, 1)); x[argmax(row_min)] = 1.0
        y = zeros(size(A, 2)); y[argmin(col_max)] = 1.0
        V = 0.5 * (maximum(row_min) + minimum(col_max))
        fallback = true
        @warn "matrix game certificate failed; falling back to pure maximin/minimax" violation=cert.violation
    end

    (; x, y, V, V_lo, V_hi,
       gap            = abs(V_hi - V_lo),
       bilinear       = game_cost(x, y, A),
       status1        = sol1.status,
       status2        = sol2.status,
       cert_violation = cert.violation,
       fallback)
end

export solve_mixed_nash

## ====================================================================

"""
    solve_mixed_security_strategy(M)

Security strategy of the row player who MINIMIZES `M`, i.e. `argmin_x max_j (x'M)_j`,
together with the value of that game.

Returns `(; x, v, status, obj)`. Note `v` (lowercase) is the *game value*, while
`solve_simplex_lp` returns `obj`, the raw LP objective `1'z` — a different
quantity. Those two used to be `v` and `V`, one character apart, which is how the
value came to be discarded at the call site.
"""
function solve_mixed_security_strategy(M_in)

    # Shift the matrix entrywise positive so the standard LP transformation applies.
    # (Previously this was guarded by `if min_value <= 0` with `min_value` seeded at
    # 0 and only ever decreasing -- so the branch always fired, and had it not,
    # `M` and `c` would have been undefined.)
    c = 1 - minimum(M_in)
    M = M_in .+ c

    ans     = solve_simplex_lp(M)
    x_tilde = ans.x

    # Transform back to the probability simplex. Normalizing by `sum(x_tilde)`
    # rather than by the solver's reported objective makes `sum(x) == 1` hold to
    # machine precision regardless of solver slop, and keeps the value consistent
    # with the strategy actually returned.
    s      = sum(x_tilde)
    x_star = s > 1e-12 ? x_tilde ./ s : fill(inv(length(x_tilde)), length(x_tilde))
    V_star = (s > 1e-12 ? 1 / s : 0.0) - c

    (; x=x_star, v=V_star, status=ans.status, obj=ans.obj)
end

export solve_mixed_security_strategy

## ====================================================================

"""
    solve_simplex_lp(A)

Solves `max 1'z s.t. A'z <= 1, z >= 0` — the standard LP of the row player who
minimizes `A`. Returns `(; x, obj, status)` and never throws.

Two changes from the original:

  * **The `z[i] >= 1e-4` floor is gone**, replaced by plain nonnegativity. That
    floor bounded `z`, not the probability `x = z/sum(z)`, so the probability
    floor it induced varied per game step -- it was not the "probability floor"
    the paper describes. Worse, it structurally forbids sparse-support
    equilibria, which most 6x6 games have, so the LP was not computing the
    minimax that the safety claim rests on.

  * **The solver is Ipopt, not OSQP.** OSQP is a first-order ADMM QP solver being
    applied to a pure LP. Measured over 300 random 6x6 games, by the saddle-point
    residual of `nash_certificate`:

        OSQP, defaults          median 1.7e-3    100% above 1e-8
        OSQP, polish + 1e-9     median 2.6e-16    24% above 1e-8   (bimodal)
        OSQP, polish + 1e-12    median 2.6e-16    25% above 1e-8   (bimodal)
        Ipopt, tol = 1e-12      median 2.8e-12     0% above 1e-8
        Ipopt, tol = 1e-14      median 2.7e-14     0% above 1e-8   <- chosen

    OSQP with polish is exact when polish succeeds and ~5e-3 when it does not,
    and tightening tolerances does not change that split. Ipopt is uniformly
    accurate, and at tol=1e-14 it is also the fastest of the three (3.3 ms vs
    9.4 ms). Ipopt is already a dependency, so this costs no Manifest change.

    This measurement is what Reviewer 7 #1 asked for, and it is why the CDC
    numbers were computed on LPs solved to roughly three decimal places.
"""
function solve_simplex_lp(A; tol=1e-14)

    model = JuMP.Model()
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_silent(model)
    JuMP.set_optimizer_attribute(model, "print_level", 0)
    JuMP.set_optimizer_attribute(model, "tol", tol)
    JuMP.set_optimizer_attribute(model, "constr_viol_tol", tol)
    JuMP.set_optimizer_attribute(model, "dual_inf_tol", tol)
    JuMP.set_optimizer_attribute(model, "compl_inf_tol", tol)
    JuMP.set_optimizer_attribute(model, "honor_original_bounds", "yes")

    r, p = size(A)

    @variable(model, z[1:r] >= 0)
    @objective(model, Max, ones(r)' * z)
    @constraint(model, c1, ones(p) >= A' * z)

    JuMP.optimize!(model)

    status = JuMP.termination_status(model)
    (; x=JuMP.value.(z), obj=JuMP.objective_value(model), status)
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
        t_horizon = params.t_horizon
        n_seg_horizon = params.n_seg_horizon
        mu = params.mu

        for jj in eachindex(vertices)

            rv_f = [vertices[jj]; v_f]
            sol = min_Δv_dist_solve(rv_0, rv_f, t_horizon, n_seg_horizon, mu;
                                    Δv_max=params.Δv_max, w_miss=params.w_miss)
            Δv_sol = sol.Δv_sol
            t, rv_hist = prop_kepler_tof_Nseg(rv_0, Δv_sol, n_seg_horizon, t_horizon / n_seg_horizon, mu)

            # save hist
            push!(p.X, rv_hist)
            push!(p.U, Δv_sol)
            push!(p.t, t)

            # Per-solve diagnostics. `miss` is the A3 reachability metric — the
            # distance between where the candidate actually ends up and the hexagon
            # vertex it was aimed at — recorded here so it is monitored continuously
            # rather than only by an offline probe. `converged` / `iters` are what
            # make the solver's convergence rate over a whole sweep reportable
            # (Reviewer 7 #1); `escalated` records whether the ΔV cap actually bound
            # and forced the augmented-Lagrangian path.
            push!(p.solve_info.traj, (
                vertex      = jj,
                t_solve     = sol.t,
                miss        = norm(rv_hist[end, 1:3] - vertices[jj]),
                max_Δv      = sol.max_Δv,
                sum_Δv      = sum(norm.(eachrow(Δv_sol))),
                converged   = sol.converged,
                g_converged = sol.g_converged,
                iters       = sol.iters,
                escalated   = sol.escalated,
            ))

        end

    end

    return players
end

export compute_players_XU


## ====================================================================

"""
    stage_cost(x1, x2, u1, u2; λ1, λ2)

P1's (the evader's) stage payoff. **P1 MAXIMIZES this**; P2's cost is its negation.

    λ1 * sqrt(‖r1 - r2‖ + 0.1)  +  λ2 * (‖u_P‖ - ‖u_E‖)

with `u1 = u_E` (evader control) and `u2 = u_P` (pursuer control). Separation up raises P1's
payoff; the evader burning fuel lowers it; the pursuer burning fuel raises it.

**Defect D4.** This previously computed `0.1 * norm(u1 - u2)` — the norm of the vector
*difference* of the two controls. That is nonnegative, so it was added to **both** players'
costs rather than transferred between them, breaking the zero-sum interpretation of the fuel
term; and it is minimized when both players thrust *identically*, which is not a fuel incentive
at all.

**Sign trap.** The "intended" form left in a comment was `0.1*(norm(u1) - norm(u2))`, which is
*inverted*: under the corrected D1 orientation this is P1's payoff and P1 maximizes it, so the
paper's Eq. (9) `λ₂(‖u_P‖ − ‖u_E‖)` requires `norm(u2) - norm(u1)`. Restoring the commented
line verbatim would have survived the A1 fix and quietly corrupted Gate B. Pinned by the
"stage cost orientation" testset.

Note the two terms are dimensionally incommensurable — `sqrt(km)` against `km/s` — so λ1 and λ2
are arbitrary scale factors, not physical weights. Say so in the paper.
"""
function stage_cost(x1, x2, u1, u2; λ1=1.0, λ2=0.1)

    dist = norm(x1[1:3] - x2[1:3])
    cost = λ1 * sqrt(dist + 0.1) + λ2 * (norm(u2) - norm(u1))

    # capture bonus, disabled: there is no capture/termination condition anywhere in
    # the codebase yet (see Gate B, B3)
    # if dist < capture_threshold
    #     cost -= 50.0
    # end

    return cost
end

export stage_cost


## ==================================================================== 
# zero-sum game 

function compute_cost_matrices(players; λ1=1.0, λ2=0.1)

    # loop through time corresponding with control inputs 
    U_idx = eachindex(players[1].U[1][:, 1])

    # start with vertex 1 for players 1 and 2 
    i_vert = 1
    j_vert = 1

    # arity is derived, not hardcoded: `n_vertices` below already computes it
    # correctly, so a 12- or 24-vertex polygon (C4) does not BoundsError here
    n_vertices_alloc = length(players[begin].t)
    player1_cost_matrix = zeros(n_vertices_alloc, n_vertices_alloc)
    player2_cost_matrix = zeros(n_vertices_alloc, n_vertices_alloc)

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
                cost1 = stage_cost(x1, x2, u1, u2; λ1, λ2)
                cost2 = -cost1
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

        # compute stage cost at each time step for each trajectory.
        # `get` guards against games saved before λ1/λ2 were added to params.
        λ1, λ2 = get(params, :λ1, 1.0), get(params, :λ2, 0.1)
        stage_cost_game = [stage_cost(x1[ii, :], x2[ii, :], u1[ii, :], u2[ii, :]; λ1, λ2) for ii in 1:size(u1, 1)]
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
    t_lp = @elapsed sol = solve_mixed_nash(players[1].cost)

    # `sol.x` and `sol.y` are already projected onto the probability simplex by
    # solve_mixed_nash. That normalization is load-bearing well beyond sampling:
    # update_strategy_belief! uses weights[chosen] as a likelihood, and
    # predict_opponent_vertices copies weights straight into a probability
    # vector. Both are meaningless on raw LP output.
    mixing_weights = (; sol.x, sol.y)
    players[1].weights = mixing_weights[1]
    players[2].weights = mixing_weights[2]

    # Retain the LP diagnostics: the Nash value used to be computed and thrown
    # away, so there was no way to check maximin == minimax, and no record of
    # whether a solve had degraded.
    lp_record = (
        V              = sol.V,
        V_lo           = sol.V_lo,
        V_hi           = sol.V_hi,
        gap            = sol.gap,
        cert_violation = sol.cert_violation,
        fallback       = sol.fallback,
        status1        = sol.status1,
        status2        = sol.status2,
        t_lp           = t_lp,
    )
    for p in players
        p.solve_info = (; traj=p.solve_info.traj, lp=lp_record)
    end

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
        # D5: two defects here, both fixed. (1) `predict_opponent_vertices` returns a
        # distribution over P1's vertices, and cost is indexed [P1_vertex, P2_vertex],
        # so the transpose is required -- without it the contraction ran over the wrong
        # axis and returned a vector indexed by P1's vertices, which was then assigned
        # to chosen[2]. Both dimensions are 6, so it silently produced a plausible
        # index instead of erroring. (2) `argmin` was wrong: players[2].cost = -A, so
        # maximizing it is what minimizes separation. Now mirrors FP_greedy exactly.
        predicted_v_probs = predict_opponent_vertices(players[2], players[1])
        expected_costs = players[2].cost' * predicted_v_probs
        chosen[2] = argmax(expected_costs)

    elseif p2_strategy == "Meta_mixed"
        # D5: transpose added, and the weight formula flipped to match FP_mixed --
        # it previously used `maximum(...) .- ec`, decreasing in expected cost, i.e.
        # the opposite of the correct branch directly above it.
        predicted_v_probs = predict_opponent_vertices(players[2], players[1])
        meta_expected_costs = players[2].cost' * predicted_v_probs
        w = meta_expected_costs .- minimum(meta_expected_costs) .+ 1e-6
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
        players[i].t_chosen = players[i].t[idx][1:params.n_seg_per_game_step+1]
        players[i].rv_chosen = players[i].X[idx][1:params.n_seg_per_game_step+1, :]
        players[i].U_chosen = players[i].U[idx][1:params.n_seg_per_game_step, :]
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
    players = compute_cost_matrices(players; λ1=params.λ1, λ2=params.λ2)

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

