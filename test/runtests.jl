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
    # Currently returns squares; flip to @test when obj_fns.jl:166 is restored.
    @test_broken sum_norm_Δv(vec([3.0 4.0 0.0; 0.0 0.0 0.0]), 2) ≈ 5.0
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

    # Differential fuel: the evader burning should LOWER P1's payoff and the
    # pursuer burning should RAISE it. The current `norm(u1-u2)` form is
    # nonnegative and so raises the payoff in both cases. Flip to @test after A4.
    @test_broken stage_cost(x1, x2, ub, u0) < stage_cost(x1, x2, u0, u0)
    @test stage_cost(x1, x2, u0, ub) > stage_cost(x1, x2, u0, u0)
end

# ---------------------------------------------------------------------------
# C3 down-payment: fictitious play on a FROZEN matrix must reach the LP solution.
# Robinson (1951) guarantees this only for a fixed payoff matrix. Fails under the
# pre-A1 orientation, passes after -- so it is also an integration test for A1.
# ---------------------------------------------------------------------------

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
