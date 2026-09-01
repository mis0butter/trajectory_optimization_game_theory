# =============================================================================
# Gate A verification figures.
#
#   julia -t 32 --project=. test/figures_gate_a.jl
#
# Emits vector PDFs to figures/gate_a/ and prints every number it plots, so the
# figures can be cross-checked against README_plan.md notes [a]-[i]. A figure
# that silently disagrees with the prose is worse than no figure.
#
# Reads the smoke-sweep data in test/results/n10/ (25 matchups x 5 games x 10
# steps, committed in 19fb075). F3, F4 and F7 additionally re-run solvers.
# =============================================================================

using trajectory_optimization_game_theory
# CairoMakie MUST be loaded after the package: src/ pulls in GLMakie
# unconditionally, and whichever backend activates last owns `save`.
using CairoMakie
CairoMakie.activate!(type = "pdf")

using LinearAlgebra, Random, Statistics, Printf
using JuMP, OSQP, Ipopt   # F3 re-solves the saved games with each LP backend
using Optim, ForwardDiff  # F4 compares inner solvers on the trajectory subproblem
using JLD2: jldsave, load   # F7 caches its game runs so the figure can be re-rendered free

const STRATS  = ["greedy", "mixed", "random", "FP_greedy", "FP_mixed"]
const NGAMES  = 5
const NSTEPS  = 10
const ROOT    = pkgdir(trajectory_optimization_game_theory)
const OUTDIR  = joinpath(ROOT, "figures", "gate_a")
const R_POLY  = 6.378                      # km, params.R_polygon

mkpath(OUTDIR)

# ---------------------------------------------------------------- style
# Sized for a two-column paper: ~252 pt single column, ~516 pt double.
const COL1, COL2 = 252, 516

gate_a_theme() = Theme(
    fontsize = 9,
    figure_padding = 6,
    Axis = (
        xgridvisible = false, ygridvisible = true,
        ygridcolor = (:black, 0.08), ygridwidth = 0.6,
        spinewidth = 0.8, xtickwidth = 0.8, ytickwidth = 0.8,
        titlesize = 10, titlealign = :left, titlefont = :bold,
        xlabelsize = 9, ylabelsize = 9,
    ),
    Legend = (framevisible = false, labelsize = 8, patchsize = (10, 10)),
)
set_theme!(gate_a_theme())

const PAL = (
    ok   = "#2166ac",   # current / correct
    bad  = "#b2182b",   # pre-fix / failing
    warn = "#e08214",
    mute = "#7f7f7f",
    alt  = "#4d9221",
)

"""
Save at TRUE SIZE plus a PNG for viewing.

`pt_per_unit = 1` matters: CairoMakie defaults to 0.75, which silently shrinks a
`size = (516, 215)` figure to a 387x161 pt page and renders 9 pt text at ~6.75 pt.
With 1.0 the page is exactly the figure in points — 516 pt = 7.16 in, the usual
double-column width — and font sizes mean what they say.
"""
function savefig(fig, name)
    png = joinpath(OUTDIR, name * ".png")
    save(png, fig; px_per_unit = 4)      # ~288 dpi raster
    # For the paper, re-enable the vector copy — pt_per_unit = 1 keeps true size:
    #   save(joinpath(OUTDIR, name * ".pdf"), fig; pt_per_unit = 1)
    println("  wrote $(basename(png))")
    png
end

# ---------------------------------------------------------------- data
"""
Walk the 25 saved matchups once and flatten the per-solve diagnostics.

Returns `(traj, lp)`. **The LP record is deduplicated**: `compute_mixing_weights!`
writes an identical copy into *both* players, so a naive flatten counts 2,500
where only 1,250 are distinct.
"""
function load_all()
    traj = NamedTuple[]
    lp   = NamedTuple[]
    cost = NamedTuple[]        # one row per step: the 6x6 payoff matrix + weights
    for p1 in STRATS, p2 in STRATS
        gv = load_games_vec(NGAMES, NSTEPS, p1, p2)
        for (gi, g) in enumerate(gv), k in eachindex(g.p1_state)
            for (pi, st) in enumerate((g.p1_state[k], g.p2_state[k]))
                for r in st.solve_info.traj
                    push!(traj, (; matchup = "$p1|$p2", p1, p2, game = gi, step = k,
                                   player = pi, r...))
                end
            end
            # LP + payoff matrix: read from p1_state only (p2 holds an identical copy)
            st1, st2 = g.p1_state[k], g.p2_state[k]
            if st1.solve_info.lp !== nothing
                push!(lp, (; matchup = "$p1|$p2", p1, p2, game = gi, step = k,
                             st1.solve_info.lp...))
            end
            if !isempty(st1.cost) && !isempty(st1.weights) && !isempty(st2.weights)
                push!(cost, (; matchup = "$p1|$p2", p1, p2, game = gi, step = k,
                               A = st1.cost, x = st1.weights, y = st2.weights))
            end
        end
    end
    (traj, lp, cost)
end

# ---------------------------------------------------------------- helpers
"ECDF points for a positive-valued sample, ready for a log-x stairs plot."
function ecdf_pts(v)
    s = sort(filter(x -> isfinite(x) && x > 0, v))
    isempty(s) && return (Float64[], Float64[])
    (s, range(1 / length(s), 1.0; length = length(s)) |> collect)
end

"Print an expected-vs-measured check line."
function check(label, measured, expected; rtol = 0.15, fmt = "%.3e")
    okmark = (isfinite(expected) && isfinite(measured) &&
              isapprox(measured, expected; rtol)) ? "ok " : "**DIFFERS**"
    @printf("  %-42s %12s   expected %-12s %s\n", label,
            Printf.format(Printf.Format(fmt), measured),
            Printf.format(Printf.Format(fmt), expected), okmark)
end

println("="^78)
println("Loading smoke-sweep data from test/results/n$(NSTEPS)/ ...")
traj, lp, cost = load_all()
@printf("  %d trajectory records, %d unique LP records, %d cost matrices\n",
        length(traj), length(lp), length(cost))

# ---------------------------------------------------------------- verification
println("\n" * "="^78)
println("VERIFICATION — measured here vs README_plan.md notes [a],[b],[f],[i]")
println("="^78)

miss   = [r.miss for r in traj]
# some solves hit the vertex exactly; log axes need a positive floor
const MISS_FLOOR = minimum(filter(>(0), miss))
maxdv  = [r.max_Δv for r in traj]
conv   = count(r -> r.converged, traj) / length(traj)
esc    = count(r -> r.escalated, traj) / length(traj)
certv  = [r.cert_violation for r in lp]
gaps   = [r.gap for r in lp]
nfall  = count(r -> r.fallback, lp)

check("terminal miss, median [km]",        median(miss),   7.6e-12)
check("terminal miss, max [km]",           maximum(miss),  8.2e-05)
@printf("  %-42s %12s   (%.1f%% of R/2)\n", "terminal miss, p97 [km]",
        Printf.format(Printf.Format("%.3e"), quantile(miss, 0.97)),
        100 * quantile(miss, 0.97) / (R_POLY/2))
@printf("  %-42s %12s   (%.1f%% of R/2)\n", "terminal miss, p99.7 (3-sigma) [km]",
        Printf.format(Printf.Format("%.3e"), quantile(miss, 0.9973)),
        100 * quantile(miss, 0.9973) / (R_POLY/2))
check("trajectory solves converged [%]",   100conv,        99.53;  fmt = "%.2f")
check("escalated to aug-L [%]",            100esc,         0.0;    fmt = "%.2f")
check("LP certificate violation, max",     maximum(certv), 1.1e-12)
check("LP fallbacks [count]",              nfall,          0.0;    fmt = "%.0f")
check("maximin-minimax gap, max",          maximum(gaps),  7.5e-13)

# NOTE: note [b] records 0.0221 km/s, measured on Day 0 under Nelder-Mead with the
# SQUARED-norm fuel term. D11 (true ΔV) and A6 (BFGS) both changed the trajectories,
# and the current value is 0.0525. The conclusion is unchanged — nothing binds — but
# the margin to the paper's stated 0.1 km/s cap is now ~2x, not ~5x.
check("max per-segment |dv| [km/s]",       maximum(maxdv), 0.0525; fmt = "%.4f")

# mixed/mixed NashConv of the strategies actually played
mm = [exploitability(c.A, c.x, c.y).nashconv for c in cost if c.p1 == "mixed" && c.p2 == "mixed"]
check("mixed/mixed NashConv, median",      median(mm),     1.6e-14)

# note [f]: matrix structure, normalized by |mean(A)|
function struct_metrics(A)
    sc = abs(mean(A)); q = fill(1 / size(A, 2), size(A, 2))
    aq = A * q; srt = sort(aq, rev = true)
    (rng = (maximum(A) - minimum(A)) / sc,
     lev = (maximum(aq) - minimum(aq)) / sc,
     gap = (srt[1] - srt[2]) / sc)
end
sm = [struct_metrics(c.A) for c in cost]
# NOTE: note [f] records 46.0 / 15.4 / 1.70, measured on 240 matrices BEFORE D4 changed
# the stage cost to differential fuel. These 1,250 matrices carry the corrected cost, so
# they supersede it. The qualitative conclusion holds — the game is not degenerate, and
# the decision margin is thin — but the margin is 2.37%, not 1.7%.
check("matrix range, median [%]",          100median(r.rng for r in sm), 59.8; fmt = "%.1f")
check("P1 leverage, median [%]",           100median(r.lev for r in sm), 18.1; fmt = "%.1f")
check("best-vs-2nd-best gap, median [%]",  100median(r.gap for r in sm), 2.37; fmt = "%.2f")
println("="^78)

# =============================================================================
# F1 — A3 reachability: the existential check
# =============================================================================
println("\nF1  reachability ...")
let
    fig = Figure(size = (COL2, 215))

    # Histogram (where the mass sits) BEHIND the ECDF (what fraction is below x).
    # Two axes in the same cell: the histogram axis is created first so it renders
    # underneath, and the ECDF axis on top is given a transparent background.
    x, y = ecdf_pts(miss)
    mpos = filter(>(0), miss)
    xlo, xhi = minimum(x) / 3, R_POLY * 2

    axh = Axis(fig[1, 1], xscale = log10, yaxisposition = :right,
               ylabel = "solves per bin", ygridvisible = false, xgridvisible = false)
    hidexdecorations!(axh); hidespines!(axh, :l, :t, :b)
    hist!(axh, mpos;
          bins = 10 .^ range(log10(minimum(mpos)), log10(maximum(mpos)); length = 42),
          color = (PAL.ok, 0.20), strokewidth = 0)
    xlims!(axh, xlo, xhi)

    ax = Axis(fig[1, 1], xscale = log10, title = "(a)  terminal miss, all 15,000 solves",
              xlabel = "‖endpoint − target vertex‖  [km]", ylabel = "cumulative fraction",
              backgroundcolor = :transparent)
    stairs!(ax, x, y, color = PAL.ok, linewidth = 1.8, step = :post)
    vlines!(ax, [R_POLY / 2], color = PAL.bad, linestyle = :dash, linewidth = 1.2)
    text!(ax, R_POLY / 2, 0.40; text = "identifiability\nthreshold R/2 →",
          align = (:right, :center), fontsize = 7, color = PAL.bad)

    # Coverage threshold. Reported as a QUANTILE, not mean ± 3sd: this distribution
    # is heavy-tailed and spans ~8 decades, so its standard deviation is set almost
    # entirely by the tail and "3 sigma" would not describe any real coverage.
    q997 = quantile(miss, 0.9973)          # the 3-sigma-equivalent coverage
    vlines!(ax, [q997], color = PAL.mute, linestyle = :dashdot, linewidth = 1.2)
    text!(ax, q997, 0.16; text = @sprintf(" 99.7%% of solves\n below %.1e km", q997),
          align = (:left, :center), fontsize = 7, color = :black)
    text!(ax, 0.03, 0.97;
          text = @sprintf("median %.1e km\np95      %.1e\nmax      %.1e", median(miss),
                          quantile(miss, 0.95), maximum(miss)),
          align = (:left, :top), space = :relative, fontsize = 7)
    xlims!(ax, xlo, xhi); ylims!(ax, -0.03, 1.05)
    linkxaxes!(ax, axh)

    ax2 = Axis(fig[1, 2], yscale = log10, title = "(b)  by game step",
               xlabel = "game step", ylabel = "terminal miss  [km]", xticks = 1:NSTEPS)
    steps = [r.step for r in traj]
    boxplot!(ax2, steps, max.(miss, MISS_FLOOR); width = 0.6, color = (PAL.ok, 0.45),
             strokecolor = PAL.ok, strokewidth = 0.6, markersize = 1.2,
             mediancolor = :black, medianlinewidth = 1.2)
    hlines!(ax2, [R_POLY / 2], color = PAL.bad, linestyle = :dash, linewidth = 1.2)
    text!(ax2, 0.98, 0.88; text = "R/2 — no late-game degradation",
          align = (:right, :top), space = :relative, fontsize = 7, color = PAL.bad)
    ylims!(ax2, MISS_FLOOR / 5, R_POLY * 3)

    colsize!(fig.layout, 1, Relative(0.45))
    savefig(fig, "fig1_reachability")
end

# =============================================================================
# F2 — D1: the LP before and after, on identical saved payoff matrices
# =============================================================================
println("F2  LP correctness (recomputing the pre-fix convention) ...")
let
    sms = solve_mixed_security_strategy
    nc_old = Float64[]; nc_new = Float64[]
    for c in cost
        A = c.A
        xo = normalize_simplex(sms(A).x)        # the ORIGINAL pairing
        yo = normalize_simplex(sms(-A').x)
        push!(nc_old, exploitability(A, xo, yo).nashconv)
        push!(nc_new, exploitability(A, c.x, c.y).nashconv)   # as actually played
    end
    @printf("  NashConv  pre-fix %.3e  |  current %.3e  |  ratio %.2e x\n",
            median(nc_old), median(nc_new), median(nc_old) / median(nc_new))
    @printf("  pre-fix solves that were NOT an equilibrium (>1e-8): %d / %d (%.1f%%)\n",
            count(>(1e-8), nc_old), length(nc_old), 100count(>(1e-8), nc_old)/length(nc_old))

    fig = Figure(size = (COL2, 215))
    ax = Axis(fig[1, 1], xscale = log10,
              title = "(a)  NashConv of the mixed strategies",
              xlabel = "NashConv   (0 ⟺ Nash equilibrium)", ylabel = "cumulative fraction")
    for (v, col, lab) in ((nc_old, PAL.bad, "pre-fix:  sms(A), sms(−Aᵀ)"),
                          (nc_new, PAL.ok,  "current:  sms(−A), sms(Aᵀ)"))
        xx, yy = ecdf_pts(max.(v, 1e-18))
        stairs!(ax, xx, yy, color = col, linewidth = 1.6, step = :post, label = lab)
    end
    vlines!(ax, [1e-8], color = PAL.mute, linestyle = :dot, linewidth = 1.0)
    text!(ax, 1e-8, 0.02; text = " 1e−8", align = (:left, :bottom), fontsize = 7,
          color = PAL.mute)
    axislegend(ax, position = :rb)
    text!(ax, 0.30, 0.78;
          text = @sprintf("medians differ\nby %.0e ×\n\n%.0f%% of pre-fix solves\nwere not equilibria",
                          median(nc_old)/median(nc_new), 100count(>(1e-8), nc_old)/length(nc_old)),
          align = (:left, :top), space = :relative, fontsize = 7)
    ylims!(ax, -0.03, 1.05)

    ax2 = Axis(fig[1, 2], title = "(b)  maximin = minimax", aspect = 1,
               xlabel = "maximin  V_lo", ylabel = "minimax  V_hi")
    lo = [r.V_lo for r in lp]; hi = [r.V_hi for r in lp]
    ablines!(ax2, 0, 1, color = PAL.mute, linestyle = :dash, linewidth = 1.0)
    scatter!(ax2, lo, hi, color = (PAL.ok, 0.35), markersize = 3, strokewidth = 0)
    text!(ax2, 0.04, 0.96; text = @sprintf("max |V_hi−V_lo|\n= %.1e   (n=%d)",
                                           maximum(gaps), length(lp)),
          align = (:left, :top), space = :relative, fontsize = 7)

    colsize!(fig.layout, 1, Relative(0.60))
    savefig(fig, "fig2_lp_correctness")
end

# =============================================================================
# F5 — D6: the ΔV cap is inert at every plausible value
# =============================================================================
println("F5  ΔV cap ...")
let
    fig = Figure(size = (COL1, 175))
    ax = Axis(fig[1, 1], xscale = log10, title = "max per-segment ‖Δv‖",
              xlabel = "‖Δv‖  [km/s]", ylabel = "solves")
    dvpos = filter(>(0), maxdv)
    hist!(ax, dvpos; bins = 10 .^ range(log10(minimum(dvpos)), log10(maximum(dvpos)); length = 45),
          color = (PAL.ok, 0.65), strokewidth = 0)
    # stagger the label heights so they do not collide
    for (v, lab, col, h) in ((maximum(maxdv), "observed max", PAL.alt,  0.97),
                             (0.1,  "paper's cap",  PAL.warn, 0.62),
                             (2.0,  "code default", PAL.bad,  0.97))
        vlines!(ax, [v], color = col, linestyle = :dash, linewidth = 1.2)
        text!(ax, v, h; text = " " * lab, align = (:left, :top), fontsize = 7,
              color = col, space = :relative, offset = (2, 0))
    end
    xlims!(ax, minimum(dvpos) / 2, 5.0)
    savefig(fig, "fig5_dv_cap")
end

# =============================================================================
# F6 — note [f]: the game is not degenerate, but the decision margin is thin
# =============================================================================
println("F6  matrix structure ...")
let
    fig = Figure(size = (COL2, 235))

    ax = Axis(fig[1, 1], title = "(a)  matrix structure — all 1,250 steps",
              ylabel = "% of |mean A|",
              xticks = (1:3, ["range", "leverage", "best vs\n2nd"]))
    for (i, f) in enumerate((:rng, :lev, :gap))
        v = 100 .* [getfield(r, f) for r in sm]
        col = i == 3 ? PAL.bad : PAL.ok
        boxplot!(ax, fill(i, length(v)), v; width = 0.5, color = (col, 0.45),
                 strokecolor = col, strokewidth = 0.6, markersize = 1.2,
                 mediancolor = :black, medianlinewidth = 1.2)
        # medians spelled out: on a linear axis the decision-margin box is a sliver,
        # which is the point, but the number still has to be readable
        text!(ax, i, median(v); text = @sprintf("  %.1f%%", median(v)),
              align = (:left, :center), fontsize = 7,
              color = i == 3 ? PAL.bad : :black)
    end
    # headroom above the range whisker for the caption, and to the right for the
    # median labels (the 2.4% one otherwise runs off the panel)
    xlims!(ax, 0.45, 3.85); ylims!(ax, -6, 132)
    text!(ax, 0.5, 0.98;
          text = "choice matters (leverage 18%),\nbut the top two are nearly tied",
          align = (:center, :top), space = :relative, fontsize = 6.5, color = PAL.mute)

    c = cost[findfirst(c -> c.p1 == "mixed" && c.p2 == "mixed" && c.step == 5, cost)]
    ax2 = Axis(fig[1, 2], title = "(b)  one payoff matrix A[i,j]", aspect = 1,
               xlabel = "P1 vertex i", ylabel = "P2 vertex j", xticks = 1:6, yticks = 1:6)
    hm = heatmap!(ax2, 1:6, 1:6, c.A, colormap = :viridis)
    Colorbar(fig[1, 3], hm, width = 7, ticklabelsize = 7)

    ax3 = Axis(fig[1, 4], title = "(c)  its equilibrium", xlabel = "weight",
               ylabel = "vertex", yticks = 1:6, xticks = -1:0.5:1)
    xw, yw = normalize_simplex(c.x), normalize_simplex(c.y)
    barplot!(ax3, 1:6,  xw; direction = :x, color = (PAL.ok, 0.85))
    barplot!(ax3, 1:6, -yw; direction = :x, color = (PAL.bad, 0.85))
    vlines!(ax3, [0], color = :black, linewidth = 0.8)
    text!(ax3, 0.95, 0.97; text = "P1", align = (:right, :top), space = :relative,
          fontsize = 7, color = PAL.ok)
    text!(ax3, 0.05, 0.97; text = "P2", align = (:left, :top), space = :relative,
          fontsize = 7, color = PAL.bad)
    xlims!(ax3, -1.05, 1.05)
    @printf("  representative equilibrium: P1 support %d/6, P2 support %d/6 (sparse ⟺ 1e-4 floor gone)\n",
            count(>(1e-6), xw), count(>(1e-6), yw))

    colsize!(fig.layout, 1, Relative(0.28)); colsize!(fig.layout, 4, Relative(0.22))
    savefig(fig, "fig6_matrix_structure")
end

# =============================================================================
# F3 — LP backend characterization  (Reviewer 7 #1)
#
# Run over the REAL saved payoff matrices, not randn(6,6): the Gate A benchmark
# used random matrices, but the actual A_k are the defensible sample and are
# already on disk.
# =============================================================================
println("F3  LP backends on real A_k ...")
let
    sub = cost[1:4:end]                       # ~310 real matrices, enough for an ECDF
    println("  benchmarking $(length(sub)) real payoff matrices x 3 backends ...")

    function lp_with(opt, A; attrs = Dict())
        m = JuMP.Model(); JuMP.set_optimizer(m, opt); JuMP.set_silent(m)
        for (k, v) in attrs; JuMP.set_optimizer_attribute(m, k, v); end
        r, p = size(A)
        JuMP.@variable(m, z[1:r] >= 0)
        JuMP.@objective(m, Max, ones(r)' * z)
        JuMP.@constraint(m, ones(p) >= A' * z)
        JuMP.optimize!(m)
        JuMP.value.(z)
    end
    function nash_with(lp, A)
        sms(M) = (c = 1 - minimum(M); zt = lp(M .+ c); s = sum(zt);
                  (x = s > 1e-12 ? zt ./ s : fill(inv(length(zt)), length(zt)),
                   v = (s > 1e-12 ? 1/s : 0.0) - c))
        s1, s2 = sms(-A), sms(A')
        V = 0.5 * (-s1.v + s2.v)
        nash_certificate(A, s1.x, s2.x, V).violation
    end

    backends = [
        ("OSQP, defaults",        PAL.bad,  A -> lp_with(OSQP.Optimizer, A)),
        ("OSQP, polish + 1e−9",   PAL.warn, A -> lp_with(OSQP.Optimizer, A;
            attrs = Dict("polish"=>true, "eps_abs"=>1e-9, "eps_rel"=>1e-9, "max_iter"=>100_000))),
        ("Ipopt, tol 1e−14",      PAL.ok,   A -> lp_with(Ipopt.Optimizer, A;
            attrs = Dict("print_level"=>0, "tol"=>1e-14, "constr_viol_tol"=>1e-14,
                         "dual_inf_tol"=>1e-14, "compl_inf_tol"=>1e-14,
                         "honor_original_bounds"=>"yes"))),
    ]

    fig = Figure(size = (COL1 + 40, 190))
    ax = Axis(fig[1, 1], xscale = log10, title = "matrix-game LP: saddle-point residual",
              xlabel = "nash_certificate violation", ylabel = "cumulative fraction")
    for (name, col, lp_fn) in backends
        v = [nash_with(lp_fn, c.A) for c in sub]
        @printf("  %-22s median %.2e   p95 %.2e   max %.2e   above 1e-8: %5.1f%%\n",
                name, median(v), quantile(v, 0.95), maximum(v), 100count(>(1e-8), v)/length(v))
        xx, yy = ecdf_pts(max.(v, 1e-18))
        stairs!(ax, xx, yy, color = col, linewidth = 1.6, step = :post, label = name)
    end
    vlines!(ax, [1e-8], color = :black, linestyle = :dot, linewidth = 1.0)
    text!(ax, 1e-8, 0.5; text = " acceptance\n threshold", align = (:left, :center),
          fontsize = 7)
    axislegend(ax, position = :lt)
    ylims!(ax, -0.03, 1.05)
    savefig(fig, "fig3_lp_solver")
end

# =============================================================================
# F4 — trajectory solver characterization + runtime  (Reviewer 7 #1)
#
# Subproblems are sampled across ALL game steps, not just step 1. The Gate A
# comparison used step-1 subproblems only and under-estimated in-pipeline cost
# by 3.3x (0.154 s measured vs 0.505 s actual); panel (c) exposes that.
# =============================================================================
println("F4  trajectory solvers, sampled across all game steps ...")
let
    params, _, _ = init_game(MersenneTwister(1))
    N, mu, th = params.n_seg_horizon, params.mu, params.t_horizon

    # rebuild (rv_0, rv_f) subproblems from a saved game, spanning every step
    subs = NamedTuple[]
    gv = load_games_vec(NGAMES, NSTEPS, "mixed", "mixed")
    for g in gv[1:3], k in eachindex(g.p1_state)
        V = collect(polygon_vertices(g.rv_ref_E[k][end, :], params))
        for (pi, st) in enumerate((g.p1_state[k], g.p2_state[k]))
            isempty(st.rv_0_hist) && continue
            j = 1 + (k + pi) % 6                       # vary the vertex across samples
            push!(subs, (step = k, rv0 = st.rv_0_hist[1, :],
                         rvf = [V[j]; st.rv_0_hist[end, 4:6]], vtx = V[j]))
        end
    end
    println("  $(length(subs)) subproblems over steps $(minimum(s.step for s in subs))–$(maximum(s.step for s in subs))")

    BT = Optim.LineSearches.BackTracking()
    function run_optim(sb, method)
        r = min_Δv_dist_solve(sb.rv0, sb.rvf, th, N, mu;
                              Δv_max = params.Δv_max, w_miss = params.w_miss, method)
        _, rh = prop_kepler_tof_Nseg(sb.rv0, r.Δv_sol, N, th / N, mu)
        (miss = norm(rh[end, 1:3] - sb.vtx), dv = sum(norm.(eachrow(r.Δv_sol))),
         t = r.t, iters = r.iters, conv = r.converged, step = sb.step)
    end
    function run_ipopt(sb)                       # scaled JuMP/@operator formulation
        tof_N, dv0 = lambert_init_guess(sb.rv0, sb.rvf, th, N, mu, "pro")
        x0 = vec(dv0); sc = max(maximum(abs, x0), 1e-6)
        raw(x) = sum_norm_Δv(x, N) +
                 params.w_miss * miss_distance_prop_kepler_Nseg(sb.rv0, x, N, sb.rvf, tof_N, mu)
        f0 = max(abs(raw(x0)), 1e-12)
        obj(z) = (v = raw(sc .* z) / f0; isfinite(v) ? v : 1e6)
        z0 = x0 ./ sc
        m = JuMP.Model(Ipopt.Optimizer); JuMP.set_silent(m)
        for (k, v) in ("print_level"=>0, "hessian_approximation"=>"limited-memory",
                       "max_iter"=>200, "bound_relax_factor"=>0.0)
            JuMP.set_optimizer_attribute(m, k, v)
        end
        JuMP.@variable(m, -20.0 <= z[k=1:3N] <= 20.0, start = z0[k])
        cfg = ForwardDiff.GradientConfig(obj, z0)
        fs(zs...) = obj(collect(zs))
        gs(g, zs...) = (ForwardDiff.gradient!(g, obj, collect(zs), cfg); nothing)
        JuMP.@operator(m, op, 3N, fs, gs); JuMP.@objective(m, Min, op(z...))
        t = @elapsed JuMP.optimize!(m)
        dvv = reshape(sc .* JuMP.value.(z), N, 3)
        _, rh = prop_kepler_tof_Nseg(sb.rv0, dvv, N, th / N, mu)
        (miss = norm(rh[end, 1:3] - sb.vtx), dv = sum(norm.(eachrow(dvv))), t = t,
         iters = JuMP.MOI.get(m, JuMP.MOI.BarrierIterations()),
         conv = JuMP.termination_status(m) in (JuMP.MOI.LOCALLY_SOLVED,
                                               JuMP.MOI.ALMOST_LOCALLY_SOLVED),
         step = sb.step)
    end

    solvers = [("Nelder-Mead", PAL.bad,  sb -> run_optim(sb, Optim.NelderMead())),
               ("Ipopt",       PAL.warn, run_ipopt),
               ("LBFGS",       PAL.alt,  sb -> run_optim(sb, Optim.LBFGS(linesearch = BT))),
               ("BFGS",        PAL.ok,   sb -> run_optim(sb, Optim.BFGS(linesearch = BT)))]

    res = Dict{String,Vector}()
    for (name, _, f) in solvers
        res[name] = [f(sb) for sb in subs]
        v = res[name]
        @printf("  %-12s conv %5.1f%%  miss %.2e  ΔV %.5f  %.3f s  %d iters (medians)\n",
                name, 100count(r -> r.conv, v)/length(v), median(r.miss for r in v),
                median(r.dv for r in v), median(r.t for r in v), median(r.iters for r in v))
    end

    fig = Figure(size = (COL2, 215))
    ax = Axis(fig[1, 1], yscale = log10, title = "(a)  accuracy vs fuel",
              xlabel = "total ΔV  [km/s]", ylabel = "terminal miss  [km]")
    for (name, col, _) in solvers
        v = res[name]
        scatter!(ax, [r.dv for r in v], max.([r.miss for r in v], 1e-16);
                 color = (col, 0.5), markersize = 4, strokewidth = 0, label = name)
    end
    axislegend(ax, position = :rt)

    ax2 = Axis(fig[1, 2], title = "(b)  converged", ylabel = "% of solves",
               xticks = (1:4, [n for (n, _, _) in solvers]), xticklabelrotation = pi/4)
    barplot!(ax2, 1:4, [100count(r -> r.conv, res[n])/length(res[n]) for (n,_,_) in solvers];
             color = [c for (_, c, _) in solvers])
    ylims!(ax2, 0, 108)

    ax3 = Axis(fig[1, 3], yscale = log10, title = "(c)  wall time per solve",
               ylabel = "seconds", xticks = (1:4, [n for (n, _, _) in solvers]),
               xticklabelrotation = pi/4)
    for (i, (name, col, _)) in enumerate(solvers)
        t = [r.t for r in res[name]]
        boxplot!(ax3, fill(i, length(t)), t; width = 0.6, color = (col, 0.5),
                 strokecolor = col, strokewidth = 0.6, markersize = 1.2,
                 mediancolor = :black, medianlinewidth = 1.2)
    end
    colsize!(fig.layout, 1, Relative(0.50))
    savefig(fig, "fig4_traj_solver")
end

# =============================================================================
# F7 — the geometric confound  (note [f])
#
# Hexagon vertices are NOT equally valuable for evasion: an in-plane radial
# offset converts into along-track drift that grows over an orbit, while
# out-of-plane displacement merely oscillates. If FP wins partly by discovering
# that fixed bias rather than by modelling an opponent, the headline is
# confounded — and B6b (best-fixed-vertex baseline) becomes mandatory.
# =============================================================================
println("F7  geometric confound: pure-vertex strategies ...")
let
    NSEED, NST = 6, 20
    cachef = joinpath(OUTDIR, "fig7_cache.jld2")
    local out
    if isfile(cachef)
        out = load(cachef, "out")
        println("  reusing cached runs ($(length(out)) games) — delete $cachef to re-run")
    else
        jobs = [(j, s) for j in 1:6 for s in 1:NSEED]
        out  = Vector{Any}(undef, length(jobs))
        t = @elapsed Threads.@threads for i in eachindex(jobs)
            j, seed = jobs[i]
            game, params = run_game(MersenneTwister(3000 + seed), NST, j, "mixed")
            sep = [norm(game.p1_state[k].rv_chosen[end, 1:3] -
                        game.p2_state[k].rv_chosen[end, 1:3]) for k in eachindex(game.p1_state)]
            out[i] = (vertex = j, seed = seed, mean_sep = mean(sep),
                      late_sep = mean(sep[max(1, end - 7):end]))
        end
        jldsave(cachef; out)
        @printf("  %d games (%d vertices x %d seeds x %d steps) in %.1f s — cached\n",
                length(jobs), 6, NSEED, NST, t)
    end

    # polygon_vertices order: (top, topin, botin, bot, botout, topout)
    # 1,4 are pure +/- axis_3 (out-of-plane); 2,3,5,6 carry an axis_2 (in-plane) component
    VNAME = ["1 top", "2 topin", "3 botin", "4 bot", "5 botout", "6 topout"]
    INPLANE = [false, true, true, false, true, true]
    for j in 1:6
        v = [o.late_sep for o in out if o.vertex == j]
        @printf("  vertex %-10s %s  late separation  mean %.3f  sd %.3f\n",
                VNAME[j], INPLANE[j] ? "in-plane " : "out-plane", mean(v), std(v))
    end

    fig = Figure(size = (COL2, 225))
    ax = Axis(fig[1, 1], title = "(a)  separation held by a fixed-vertex evader",
              subtitle = "$(NSEED) seeds x $(NST)-step games per vertex;  P1 plays vertex j every step,\n" *
                         "P2 samples the Nash mix afresh each step",
              subtitlesize = 7, subtitlegap = 2,
              xlabel = "hexagon vertex", ylabel = "late-game separation  [km]",
              xticks = (1:6, VNAME), xticklabelrotation = pi/5)
    for j in 1:6
        v = [o.late_sep for o in out if o.vertex == j]
        col = INPLANE[j] ? PAL.ok : PAL.warn
        boxplot!(ax, fill(j, length(v)), v; width = 0.55, color = (col, 0.45),
                 strokecolor = col, strokewidth = 0.7, markersize = 3,
                 mediancolor = :black, medianlinewidth = 1.2)
        scatter!(ax, fill(j, length(v)) .+ 0.22, v; color = (col, 0.7), markersize = 3)
    end
    text!(ax, 0.98, 0.97; text = "blue = in-plane (axis 2)\norange = out-of-plane (axis 3)",
          align = (:right, :top), space = :relative, fontsize = 7)
    # n=6 is still small: report the standard error so the spread is not over-read
    allv = [o.late_sep for o in out]
    @printf("  spread of means %.2f km on a grand mean of %.2f (SE per vertex ~%.2f)\n",
            maximum(mean(o.late_sep for o in out if o.vertex == j) for j in 1:6) -
            minimum(mean(o.late_sep for o in out if o.vertex == j) for j in 1:6),
            mean(allv), std(allv) / sqrt(NSEED))

    ax2 = Axis(fig[1, 2], title = "(b)  per-seed ranking — no stable order",
               subtitle = "each line is one vertex; mean over the last 8 steps",
               subtitlesize = 7, subtitlegap = 2,
               xlabel = "seed", ylabel = "rank  (1 = best)", yticks = 1:6, xticks = 1:NSEED)
    for j in 1:6
        rk = Int[]
        for s in 1:NSEED
            vals = [(o.vertex, o.late_sep) for o in out if o.seed == s]
            order = sortperm([v for (_, v) in vals], rev = true)
            push!(rk, findfirst(==(j), [vals[o][1] for o in order]))
        end
        col = INPLANE[j] ? PAL.ok : PAL.warn
        lines!(ax2, 1:NSEED, rk, color = (col, 0.85), linewidth = 1.4)
        scatter!(ax2, 1:NSEED, rk, color = col, markersize = 5)
        text!(ax2, NSEED + 0.12, rk[end]; text = " $j", align = (:left, :center),
              fontsize = 7, color = col)
    end
    ax2.yreversed = true
    xlims!(ax2, 0.6, NSEED + 0.8)

    colsize!(fig.layout, 1, Relative(0.55))
    savefig(fig, "fig7_vertex_geometry")
end

println("\n" * "="^78)
println("All figures written to $OUTDIR")
println("="^78)
