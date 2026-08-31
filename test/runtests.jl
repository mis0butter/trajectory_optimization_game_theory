using trajectory_optimization_game_theory
using Test
using LinearAlgebra, Random, Statistics
using ForwardDiff, FiniteDifferences
using JuMP

# ============================================================================
# Convention under test, everywhere below:
#   A is P1's payoff matrix.  P1 (the evader) MAXIMIZES A.
#   P2 (the pursuer) MINIMIZES A.
# ============================================================================

@testset "trajectory_optimization_game_theory" begin

# ---------------------------------------------------------------------------
# A1 / D1 — the matrix game solver
# ---------------------------------------------------------------------------

@testset "von Neumann: maximin == minimax" begin
    rng = MersenneTwister(0)
    for _ in 1:200
        r, p = rand(rng, 2:8), rand(rng, 2:8)
        A = 10 .* randn(rng, r, p) .+ 3          # negative entries exercise the shift
        s = solve_mixed_nash(A)
        lo = -solve_mixed_security_strategy(-A).v    # P1 maximin
        hi =  solve_mixed_security_strategy(A').v    # P2 minimax
        tol = 1e-6 * max(1, maximum(abs, A))
        @test isapprox(lo, hi;   atol = tol)
        @test isapprox(s.V, hi;  atol = tol)
        @test isapprox(s.x' * A * s.y, s.V; atol = tol)
    end
end

@testset "exploitability of the LP solution is ~0 by construction" begin
    # The single strongest check on A1: if the orientation were wrong, the
    # returned pair would not be a saddle point and nashconv would be large.
    rng = MersenneTwister(1)
    for _ in 1:300
        r, p  = rand(rng, 2:8), rand(rng, 2:8)
        A     = rand(rng, (0.01, 1.0, 100.0)) .* randn(rng, r, p)   # exercise scaling
        s     = solve_mixed_nash(A)
        e     = exploitability(A, s.x, s.y)
        scale = max(1.0, maximum(abs, A))
        @test e.gain1 >= -1e-8 * scale        # both gains nonnegative by definition
        @test e.gain2 >= -1e-8 * scale
        @test e.nashconv <= 1e-5 * scale      # <-- THE check on A1
        @test isapprox(e.v, s.V; atol = 1e-5 * scale)
    end
end

@testset "the pre-fix orientation IS exploitable (pins D1)" begin
    # Reproduces the original `sms(A)` / `sms(-A')` pairing and asserts it is
    # badly exploitable, so the bug cannot silently return.
    A = [1.0 2.0 0.0; 3.0 4.0 3.0; 0.0 1.0 -1.0]      # saddle at row 2, V = 3
    x_bad = normalize_simplex(solve_mixed_security_strategy(A).x)
    y_bad = normalize_simplex(solve_mixed_security_strategy(-A').x)
    @test exploitability(A, x_bad, y_bad).nashconv > 1.0

    # and the corrected solver is not exploitable on the same matrix
    s = solve_mixed_nash(A)
    @test exploitability(A, s.x, s.y).nashconv < 1e-6
end

@testset "known saddle point" begin
    # row minima 0, 3, -1  -> maximin 3 at row 2
    # col maxima 3, 4,  3  -> minimax 3
    A = [1.0 2.0 0.0; 3.0 4.0 3.0; 0.0 1.0 -1.0]
    s = solve_mixed_nash(A)
    @test isapprox(s.V, 3.0; atol = 1e-5)
    # P1 must play row 2. Pre-fix this returned ~[0,0,1], P1's WORST row.
    @test isapprox(s.x, [0.0, 1.0, 0.0]; atol = 1e-4)
end

@testset "matching pennies" begin
    s = solve_mixed_nash([1.0 -1.0; -1.0 1.0])
    @test isapprox(s.V, 0.0; atol = 1e-5)
    @test isapprox(s.x, [0.5, 0.5]; atol = 1e-4)
    @test isapprox(s.y, [0.5, 0.5]; atol = 1e-4)
end

@testset "sparse support survives (1e-4 floor removed)" begin
    # Row 1 dominates for the maximizer, so row 2 must get ~zero weight.
    # The old `z[i] >= 1e-4` constraint made this structurally impossible.
    s = solve_mixed_nash([5.0 5.0; 0.0 0.0])
    @test s.x[2] < 1e-6
end

@testset "affine invariance and role swap" begin
    A = randn(MersenneTwister(3), 6, 6)
    s = solve_mixed_nash(A)
    for (α, β) in ((1.0, 7.3), (2.5, 0.0), (2.5, -11.0))
        t = solve_mixed_nash(α .* A .+ β)          # value is affine-equivariant
        @test isapprox(t.V, α * s.V + β; atol = 1e-5)
    end
    u = solve_mixed_nash(-A')                      # swapping roles negates the value
    @test isapprox(u.V, -s.V; atol = 1e-5)
end

@testset "simplex validity, termination, certificate" begin
    rng = MersenneTwister(2)
    for _ in 1:200
        A = randn(rng, 6, 6)
        s = solve_mixed_nash(A)
        @test all(>=(-1e-12), s.x) && all(>=(-1e-12), s.y)
        @test isapprox(sum(s.x), 1.0; atol = 1e-10)
        @test isapprox(sum(s.y), 1.0; atol = 1e-10)
        @test s.cert_violation < 1e-8      # solver-independent; the substantive one
        @test !s.fallback                  # no solve should need the fallback here
        @test s.gap < 1e-6 * max(1, maximum(abs, A))
    end
end

@testset "nash_certificate rejects a non-equilibrium" begin
    A = [1.0 2.0 0.0; 3.0 4.0 3.0; 0.0 1.0 -1.0]
    @test nash_certificate(A, [0.0,1.0,0.0], [1.0,0.0,0.0], 3.0).violation < 1e-9
    @test nash_certificate(A, [1.0,0.0,0.0], [0.0,1.0,0.0], 3.0).violation > 1e-3
end

# ---------------------------------------------------------------------------
# A2 — exploitability module
# ---------------------------------------------------------------------------

@testset "exploitability basics" begin
    A = [1.0 2.0 0.0; 3.0 4.0 3.0; 0.0 1.0 -1.0]
    e = exploitability(A, [0.0,1.0,0.0], [1.0,0.0,0.0])     # the saddle
    @test isapprox(e.nashconv, 0.0; atol = 1e-10)
    @test isapprox(e.v, 3.0; atol = 1e-10)

    e2 = exploitability(A, [1.0,0.0,0.0], [1.0,0.0,0.0])    # P1 off-equilibrium
    @test e2.gain1 > 0                                       # P1 could improve
    @test e2.nashconv > e.nashconv

    # accepts unnormalized input (fp_belief holds raw counts)
    @test isapprox(exploitability(A, [0.0,2.0,0.0], [4.0,0.0,0.0]).nashconv, 0.0; atol = 1e-10)
end

# ---------------------------------------------------------------------------
# D9 / D10 — automatic differentiation through the trajectory objective
# ---------------------------------------------------------------------------

@testset "AD plumbing through Kepler propagation (D10)" begin
    mu = 398600.4415
    rv0 = kep2cart([6998.0, 0.01, 20*pi/180, 10*pi/180, 20*pi/180, 25*pi/180], mu)

    # cart2kep must not preallocate a Vector{Float64} and assign Duals into it
    f_kep(rv) = sum(cart2kep(rv, mu))
    g_ad = ForwardDiff.gradient(f_kep, rv0)
    g_fd = FiniteDifferences.grad(central_fdm(5, 1), f_kep, rv0)[1]
    @test all(isfinite, g_ad)
    @test isapprox(g_ad, g_fd; rtol = 1e-5)

    # and through the full objective that min_Δv_dist minimizes
    N, tof = 10, 1165.2
    rvf = rv0 .+ [5.0, 3.0, 1.0, 0.0, 0.0, 0.0]
    tof_N, dv0 = lambert_init_guess(rv0, rvf, tof, N, mu, "pro")
    obj(x) = sum_norm_Δv(x, N) + miss_distance_prop_kepler_Nseg(rv0, x, N, rvf, tof_N, mu)
    x0 = vec(dv0)
    h_ad = ForwardDiff.gradient(obj, x0)
    h_fd = FiniteDifferences.grad(central_fdm(5, 1), obj, x0)[1]
    @test all(isfinite, h_ad)
    @test isapprox(h_ad, h_fd; rtol = 1e-5)
end

@testset "flat index convention is column-major" begin
    # segment i's components live at flat indices i, N+i, 2N+i -- NOT 3(i-1)+1:3i.
    # Getting this wrong would silently constrain the wrong quantity in A6.
    N = 10
    x = vec(randn(MersenneTwister(8), N, 3))
    for i in 1:N
        @test reshape(x, N, 3)[i, :] == x[[i, N + i, 2N + i]]
    end
end

@testset "sum_norm_Δv returns norms, not squares (D3b/D11)" begin
    # [3 4 0; 0 0 0] -> ‖(3,4,0)‖ + ‖(0,0,0)‖ = 5, not 9+16 = 25.
    @test sum_norm_Δv(vec([3.0 4.0 0.0; 0.0 0.0 0.0]), 2) ≈ 5.0 atol=1e-9

    # positively homogeneous of degree 1 (a norm); the squared form was degree 2
    x = vec(randn(MersenneTwister(9), 10, 3))
    @test sum_norm_Δv(2 .* x, 10) ≈ 2 * sum_norm_Δv(x, 10) rtol=1e-9

    # the ε desingularization must keep the gradient finite at Δv = 0, where a
    # plain norm would hand ForwardDiff a NaN
    g = ForwardDiff.gradient(z -> sum_norm_Δv(z, 10), zeros(30))
    @test all(isfinite, g)
    @test sum_norm_Δv(zeros(30), 10) < 1e-9
end

# ---------------------------------------------------------------------------
# Geometry and the game pipeline
# ---------------------------------------------------------------------------

@testset "polygon geometry" begin
    params, players, game = init_game(MersenneTwister(5))
    rv = game.rv_ref_E[end][end, :]
    V  = collect(polygon_vertices(rv, params))
    R  = params.R_polygon
    a1, a2, a3 = axis_123(rv)

    # axis_123 is NOT an orthonormal frame: a2 and a3 are mutually orthogonal, but
    # a1 (the velocity direction) is tilted off a2 by roughly the eccentricity
    # (measured a1.a2 = -9.9e-3 at e = 0.01). The hexagon therefore lies in the
    # a2-a3 plane, which is the local radial/normal plane -- NOT the plane normal
    # to velocity. Those coincide only for a circular orbit; here they differ by
    # about 0.57 deg. The correct invariant is membership in span(a2, a3).
    n̂ = normalize(cross(a2, a3))
    @test abs(dot(a2, a3)) < 1e-12
    @test abs(dot(a1, a2)) > 1e-4          # documents the non-orthogonality

    @test length(V) == 6
    for i in 1:6
        @test isapprox(norm(V[i] - rv[1:3]), R; rtol = 1e-10)            # circumradius
        @test isapprox(norm(V[i] - V[mod1(i+1, 6)]), R; rtol = 1e-8)     # adjacent spacing == R
        @test abs(dot(V[i] - rv[1:3], n̂)) < 1e-10                        # in the ê2-ê3 plane
    end
end

@testset "zero-sum property and cost matrix shape" begin
    _, players, _ = init_game(MersenneTwister(4))
    @test players[2].cost == -players[1].cost
    @test size(players[1].cost) == (6, 6)
    @test all(isfinite, players[1].cost)
end

@testset "weights are genuine probability vectors" begin
    # update_strategy_belief! uses weights[chosen] as a likelihood and
    # predict_opponent_vertices copies weights into a probability vector, so
    # this must hold for those to mean anything.
    _, players, _ = init_game(MersenneTwister(4))
    for p in players
        @test isapprox(sum(p.weights), 1.0; atol = 1e-10)
        @test all(>=(-1e-12), p.weights)
    end
end

@testset "A3 reachability: candidates reach their vertices" begin
    params, players, game = init_game(MersenneTwister(6))
    V = collect(polygon_vertices(game.rv_ref_E[end][end, :], params))
    miss = [norm(players[i].X[j][end, 1:3] - V[j]) for i in 1:2, j in 1:6]
    R = params.R_polygon
    @info "A3" median = median(vec(miss)) p95 = quantile(vec(miss), 0.95) R = R
    # Measured Aug 31: median 4.5e-5 km, max 1.8e-4 km -- five orders of
    # magnitude inside the identifiability threshold.
    @test quantile(vec(miss), 0.95) < R / 2
    @test maximum(miss) < R / 100
end

@testset "solve_info is populated" begin
    _, players, _ = init_game(MersenneTwister(1))
    for p in players
        @test length(p.solve_info.traj) == 6
        @test all(r -> r.t_solve > 0, p.solve_info.traj)
        @test p.solve_info.lp !== nothing
        @test p.solve_info.lp.cert_violation < 1e-8
        @test !p.solve_info.lp.fallback
    end
end

@testset "stage cost orientation (A4/D4)" begin
    x1, x2 = [7000.0,0,0,0,7.5,0], [7010.0,0,0,0,7.5,0]
    u0, ub = [0.0,0,0], [0.1,0,0]
    far    = [7100.0,0,0,0,7.5,0]

    # separation up -> P1's payoff up. True already.
    @test stage_cost(far, x2, u0, u0) > stage_cost(x1, x2, u0, u0)

    # Differential fuel λ2*(‖u_P‖ - ‖u_E‖): the evader burning LOWERS P1's payoff,
    # the pursuer burning RAISES it. The old `norm(u1-u2)` form was nonnegative and
    # so raised it in both cases.
    @test stage_cost(x1, x2, ub, u0) < stage_cost(x1, x2, u0, u0)
    @test stage_cost(x1, x2, u0, ub) > stage_cost(x1, x2, u0, u0)

    # equal burns must cancel exactly — the defining property of a differential-fuel
    # term, and precisely what norm(u1-u2) got wrong (it penalized both players for
    # thrusting in different directions)
    @test stage_cost(x1, x2, ub, ub) ≈ stage_cost(x1, x2, u0, u0)
    @test stage_cost(x1, x2, [0.1,0,0], [0.0,0.1,0]) ≈ stage_cost(x1, x2, u0, u0)

    # λ1/λ2 are honored
    @test stage_cost(x1, x2, u0, ub; λ2=0.0) ≈ stage_cost(x1, x2, u0, u0; λ2=0.0)
    @test stage_cost(far, x2, u0, u0; λ1=2.0) > stage_cost(far, x2, u0, u0; λ1=1.0)
end

# ---------------------------------------------------------------------------
# C3 down-payment: fictitious play on a FROZEN matrix must reach the LP solution.
# Robinson (1951) guarantees this only for a fixed payoff matrix. Fails under the
# pre-A1 orientation, passes after -- so it is also an integration test for A1.
# ---------------------------------------------------------------------------

@testset "gradient-based inner solver (A6/D9)" begin
    # min_optim must actually USE the gradient now. On a quadratic it should hit the
    # exact minimum in a handful of iterations; Nelder-Mead would need dozens.
    Q = Diagonal([1.0, 10.0, 100.0, 1000.0]); b = [1.0, -2.0, 3.0, -4.0]
    f(x) = 0.5 * x' * Q * x - b' * x
    r = min_optim_info(f, zeros(4))
    @test r.converged
    @test isapprox(r.x_min, Q \ b; rtol = 1e-6)
    @test r.iters < 60                      # gradient-based, not simplex

    # shape preservation: aug_L passes a 30x1 Matrix and then takes norm(x_min - x_k)
    r2 = min_optim_info(x -> sum(abs2, x), zeros(6, 1))
    @test size(r2.x_min) == (6, 1)

    # the guard: a non-finite objective must not propagate (it trips an assertion
    # inside Optim's line search rather than just returning a bad answer)
    g(x) = (v = sum(abs2, x) - 1.0; v < 0 ? NaN : v)
    @test_nowarn min_optim_info(x -> (v = g(x); isfinite(v) ? v : 1e6), [3.0, 3.0])
end

@testset "min_Δv_dist_solve diagnostics and escalation" begin
    params, players, game = init_game(MersenneTwister(12))
    V = collect(polygon_vertices(game.rv_ref_E[end][end, :], params))
    rv0 = players[1].rv_0_hist[1, :]
    rvf = [V[1]; players[1].rv_0_hist[end, 4:6]]
    N, mu, th = params.n_seg_horizon, params.mu, params.t_horizon

    s = min_Δv_dist_solve(rv0, rvf, th, N, mu; Δv_max = params.Δv_max, w_miss = params.w_miss)
    @test s.converged
    @test size(s.Δv_sol) == (N, 3)
    @test !s.escalated                       # 2.0 km/s cap is inactive by ~90x (D6)
    @test s.max_Δv < params.Δv_max
    @test s.iters > 0 && s.t > 0

    # it must actually reach the vertex
    _, rvh = prop_kepler_tof_Nseg(rv0, s.Δv_sol, N, th / N, mu)
    @test norm(rvh[end, 1:3] - V[1]) < 1e-6

    # and a genuinely binding cap must take the augmented-Lagrangian path
    s2 = min_Δv_dist_solve(rv0, rvf, th, N, mu; Δv_max = 1e-3, w_miss = params.w_miss)
    @test s2.escalated

    # min_Δv_dist keeps its old contract: just the Δv matrix
    @test min_Δv_dist(rv0, rvf, th, N, mu; Δv_max = params.Δv_max) isa AbstractMatrix
end

@testset "Meta_* / FP_* P2 axis and sign convention (D5)" begin
    # `cost` is indexed [P1_vertex, P2_vertex]. `predict_opponent_vertices` returns a
    # distribution over the OPPONENT's vertices, so for P2 that is a distribution over
    # P1's rows -- which must contract the P1 axis, i.e. multiply by cost', not cost.
    A  = [1.0 5.0; 2.0 3.0]        # deliberately non-symmetric
    C2 = -A                         # players[2].cost
    pred = [0.0, 1.0]               # P1 certainly plays vertex 2; A[2,:] = [2,3] -> P2 wants col 1

    @test argmax(C2' * pred) == 1               # correct: transpose, and argmax on -A
    @test argmax(C2  * pred) == 2               # the old un-transposed form: wrong answer
    @test argmin(C2' * pred) == 2               # the old argmin: also wrong
    @test argmax(C2' * pred) == argmin(A' * pred)   # mirrors FP_greedy exactly

    # Meta_mixed weights must be INCREASING in P2's expected cost, like FP_mixed.
    # The old form used `maximum(ec) .- ec`, which is decreasing.
    ec = C2' * pred
    @test all(diff((ec .- minimum(ec) .+ 1e-6)[sortperm(ec)]) .>= 0)

    # end-to-end: a Meta_greedy pursuer must pick a valid column index
    _, players, _ = init_game(MersenneTwister(11), "mixed", "Meta_greedy")
    @test players[2].chosen in 1:size(players[1].cost, 2)
end

@testset "frozen-matrix FP converges to the LP equilibrium" begin
    A = randn(MersenneTwister(7), 6, 6)
    s = solve_mixed_nash(A)
    f1, f2 = ones(6), ones(6)          # belief counts, as initialized at IC.jl
    for _ in 1:20_000
        i = argmax(A  * (f1 ./ sum(f1)))     # P1 maximizes
        j = argmin(A' * (f2 ./ sum(f2)))     # P2 minimizes
        f2[i] += 1
        f1[j] += 1
    end
    p̂, q̂ = f2 ./ sum(f2), f1 ./ sum(f1)
    @test exploitability(A, p̂, q̂).nashconv < 5e-2
    @test isapprox(game_cost(p̂, q̂, A), s.V; atol = 5e-2 * max(1, maximum(abs, A)))
end

end # top-level testset
