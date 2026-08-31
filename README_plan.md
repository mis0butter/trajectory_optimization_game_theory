# AAMAS 2027 submission — handoff plan

**This document is written to be handed to an agent on a different machine.** It assumes no prior context about this repository.

---

## Status — Aug 31 2026

**Gate A is 9 of 10 items done.** Only A6 (Ipopt for the trajectory optimizer) and the
smoke sweep remain. Per-defect detail is in the §3 `Status` column and the notes after it.

**Landed (all verified value-preserving except where a fix is the point):**

- **A1 / D1 — the LP is fixed.** `sms(-A)` for P1, `sms(A')` for P2, orientation derived from
  `stage_cost` rather than assumed. Value threaded (`V`, `V_lo`, `V_hi`, `gap`); weights
  normalized centrally; the `z >= 1e-4` floor dropped; `zeros(6,6)` derived from `n_vertices`.
- **LP backend OSQP → Ipopt** (`tol=1e-14`), decided by measurement — see note [c].
- **`nash_certificate`** — solver-independent saddle-point residual, plus a never-throw
  record-and-fall-back policy.
- **A6.0 / D10 — the AD path works for the first time.** One line in `cart2kep`. Gated bitwise
  on 20,000 states and against `central_fdm(5,1)` to 1.6e-9.
- **A5 / D8** — `@exfiltrate` removed, `load_games_vec` repaired, `prop_rv_ref` aliasing bug
  fixed (`copy` on `params.kep0_ref_E`), `@elapsed` instrumentation on every trajectory solve
  and every LP.
- **D6 plumbing** — `dm`/`Δv_max` keyword-only; `params` gained `Δv_max`, `w_miss`, `λ1`, `λ2`.
- **Schema (was Gate C's C2, pulled forward)** — `player_struct` gained `solve_info` and
  `learner_state` plus an all-keyword constructor, so future fields stop requiring lockstep
  edits at three sites. Done now because there is zero saved data to invalidate.
- **`src/game_theory/exploitability.jl`** — `exploitability`, `exploitability_series`,
  `belief_exploitability`.
- **A4 / D4 — stage cost is differential fuel**, `λ1·sqrt(dist+0.1) + λ2·(‖u_P‖ − ‖u_E‖)`, with
  the sign trap in note [d] avoided. First change that deliberately moves results.
- **D11 — `sum_norm_Δv` returns true ΔV**, so fuel is an active term rather than outweighed
  ~1000:1. See note [g], including the 162%-of-Lambert result.
- **D5 — both `Meta_*` P2 branches corrected** (transpose + sign), now mirroring `FP_*` exactly.
- **`CLAUDE.md` written** (task zero). Carries a compact D1–D13 index pointing here as the source
  of truth, rather than a verbatim copy that would drift.
- **`test/runtests.jl`** — was 6 lines testing nothing; now **3092 passing, 0 failing, 0 broken**.

- **A6 / D9 — the trajectory optimizer is gradient-based.** BFGS + BackTracking, 12/12
  converged, 38% less fuel, ~10^7 better terminal miss. **Ipopt was tried and does not converge
  on this problem** — see note [h]. Zero new dependencies either way.

- **Smoke sweep passed** — see note [i]: 99.53% solver convergence, 0/2,500 LP certificate
  failures, and `mixed`/`mixed` NashConv of 1.6e-14 confirming D1 end to end.

**GATE A IS COMPLETE.** Next is Gate B, starting with B1 (real Monte Carlo / D2). Note the cost
warning in [i]: a full 5x5 re-run now projects to ~2.5 h rather than ~36 min.

**Two headline changes to the plan's premises:**

1. **D3 is dead** — see note [a]. §9 risk 2 does not fire, and A3 leaves the critical path.
2. **§9 risk 4 was wrong** — the AD path was *not* de-risked; it was broken (D10), and the
   LBFGS fallback shared the same broken dependency. Fixed, but it means A6's justification is
   now D11 + convergence statistics, not reachability.

**Corrections to §1 (this machine):**

- **`test/results/` does not exist here** and was never copied — it is gitignored, so the ~12 GB
  never left the old machine. Gate A generated its own inputs instead; this cost ~1 minute.
- **32 cores / 62 GB RAM / 1.4 TB free**, not 8 cores. Measured **56 ms** per trajectory solve
  in steady state (the 1.19 s first-call figure is JIT), **46 s** per 30-step game
  single-threaded → **~45 min for a full 5×5 sweep**, not 2 h 51 min.
- Julia is **1.12.7**; `Manifest.toml` was resolved under 1.11.1. Working so far, but any
  `Pkg.add` re-resolves it — see A6.
- The repo is on branch **`aamas`**, not `main`. `README.md` does not exist.

---

## 0. Context

`/Users/june/research/game_theory_space` is a Julia package (`trajectory_optimization_game_theory`) implementing a two-spacecraft orbital pursuit-evasion game. An evader (P1) follows a reference orbit while a pursuer (P2) tries to close distance. At each game step both players generate six candidate trajectories targeting the vertices of a hexagon around the reference orbit, form a 6×6 zero-sum matrix game, solve it by LP, and pick a vertex according to a strategy (mixed / greedy / random / fictitious-play variants).

This work was submitted to **CDC 2026 (submission #2010) and rejected.** Three reviews plus an AE summary are in `Game_Theory_Space_CDC_2026/reviews/`. They converge on one theme: the game-theoretic formulation is under-specified and internally inconsistent. Investigation showed **two of those complaints are genuine code bugs, not notation problems.**

Target: **AAMAS 2027.** Abstract **Oct 1 2026**, paper **Oct 8 2026**, 8 pages + unlimited references, ACM sigconf, **double-blind**, LaTeX mandatory. All co-authors need OpenReview accounts ~**Sep 17**. ([submission instructions](https://warwick.ac.uk/fac/sci/dcs/aamas2027/calls/instructions/))

**Division of labor: Junette Hsin owns the paper in Overleaf. This plan covers code, experiments, analysis, and figures only.** Deliver prose as markdown and figures as PDF. Do not create LaTeX in this repo.

### The agreed framing

> **Headline:** what a player gains by modeling and exploiting a specific opponent online, and what it gives up in safety by doing so.
> **Supporting lemma:** the lifting — a proposition bounding the gap to the continuous game, plus a vertex-count ablation.

The dependency runs lifting → safety: a minimax strategy over six hexagon vertices is unexploitable only against opponents confined to those six vertices, so the safety claim is worth exactly as much as the lifting is faithful.

**Positioning against the closest related work.** Peters et al., *Learning Mixed Strategies in Trajectory Games* (RSS 2022) — already cited as [3] — lifts trajectory games to finite candidate sets, computes mixed Nash per step, and plays in receding horizon (their §V-D2: 20-step horizon, replan every 9 steps, 500 updates). So lifting, mixed strategies, receding horizon, *and* a time-varying payoff matrix are all prior art. Their Table II is already a receding-horizon strategy tournament with mean ± SEM over randomized initial states.

The defensible delta is stated in their own words: *"each player is oblivious to their opponent's decision making process and solves its own version of the game."* **Peters et al. do no opponent modeling** — their adaptation is a trajectory generator trained offline by self-play, and within a match nobody infers anything about who they are playing. This paper's FP accumulates beliefs about a *specific* opponent across game steps and best-responds, with no training. That is the axis: **online opponent exploitation vs. oblivious equilibrium play.**

Do **not** claim "the game is time-varying" as novelty — it is not. Treat the time-variation instead as the thing that makes cross-step belief pooling *questionable*, and address it honestly (see C1/C3).

---

## 1. Environment setup (do this first on the new machine)

- Julia with the repo's `Project.toml` / `Manifest.toml`. `Pkg.instantiate()` from the repo root.
- **Launch Julia with threads:** `julia -t 8 --project=.`. `Threads.nthreads()` defaults to **1**, and `run_MC_games_parallel` uses `Threads.@threads` — without this every sweep runs serially.
- All result paths are **relative**; scripts must be run from the repo root.
- **Disk:** ~~`test/results/` currently holds ~12 GB.~~ **Superseded — see Status above: it does not exist on this machine.** Gate A diagnostics need only `test/results/n30/` (~1.4 GB, 25 files at ~59 MB) copied to the new machine. Budget 20–40 GB for new runs unless B2 (slim persistence) lands first.
- **Measured baseline cost (superseded — see Status):** the 5×5 sweep at 50 games × 30 steps took **2 h 51 min** wall clock on an 8-core machine, ~6 min per matchup, derived from output file mtimes.

---

## 2. Task zero — `CLAUDE.md`

**DONE Aug 31** — `CLAUDE.md` exists at the repo root. One deviation: the defect table is a
compact one-line-per-defect index pointing at §3 rather than a verbatim copy, so the two cannot
drift. (`README.md` does not exist on this machine, contrary to the note below.)

Before any code changes, write `CLAUDE.md` at the repo root. It must capture:

- **Purpose and the game:** evader/pursuer, hexagonal lifting, per-step matrix game, receding horizon.
- **Code organization:** `src/game_theory/` (matrix game solver, play loop), `src/Opt/` (trajectory optimization: augmented Lagrangian, objective/constraint functions), `src/Dyn/` (Kepler propagation), `src/Lambert/` (initial guess), `src/Utils/` (ICs, structs, plotting, MC statistics), `src/old/` (dead). Entry points: `test/run_all_MC_table.jl` produces the sweep, `test/analyze_late_game.jl` post-processes into `test/MC_results_late_game.csv` and figures.
- **Execution path for one game step**, which is the thing a new agent most needs: `prop_game_step` → `compute_states_nash` → `compute_players_XU` (12 trajectory optimizations) → `compute_cost_matrices` → `compute_mixing_weights!` (LP) → `choose_strategies!` → `update_beliefs!` → `update_chosen_trajectories!`.
- **CDC rejection and reviewer feedback**, summarized, pointing at `Game_Theory_Space_CDC_2026/reviews/`.
- **The defect table from §3 below**, verbatim — this is the highest-value content for a new agent.
- **Paper-source situation:** the LaTeX in `Game_Theory_Space_CDC_2026/` (`root.tex`, `sections/`, `refs.bib`) is the *old pre-FP draft*, not the CDC submission. The submitted version exists in this repo only as `Game_Theory_Space_ICRA_2026 (2).pdf`; its source lives in Overleaf. `test/results.md` is a third, intermediate write-up.
- **Gotchas:** threads default to 1; relative paths; `Manifest.toml` has no LP solver but OSQP; there is no real test suite (`test/runtests.jl` is 6 lines); `test/` is scripts, not tests.

---

## 3. Verified defects

D1–D9 are the original audit. **D10–D13 were found on Aug 31** while executing Gate A; the
`Status` column and the notes after the table are the current state.

| # | Defect | Location | Consequence | Status |
|---|---|---|---|---|
| D1 | **LP sign inversion.** `solve_mixed_security_strategy` returns the row player's *minimizing* security strategy. `solve_mixed_nash(A)` therefore gives the evader a distance-minimizing mixed strategy and the pursuer a distance-maximizing one — both reversed. Verified numerically on a matrix with a known saddle point: the evader is assigned ~1.0 weight on its worst row. | `matrix_game_solver.jl:55-59` | Every `mixed` and `greedy` number in the CDC tables is the equilibrium of the reversed game. Reviewer 7 #2, Reviewer 8 #1. | **FIXED** Aug 31. `sms(-A)`/`sms(A')`, orientation derived from `stage_cost`. Pinned by a test asserting the old pairing is exploitable. |
| D2 | **No Monte Carlo.** `init_game` hardcodes both orbits; `rand_IC` exists but its only call site is commented out. Only the vertex-sampling RNG varies between trials. | `IC.jl:19-26`, `IC.jl:24-25` | For deterministic matchups (greedy/greedy, greedy/FP-greedy, FP-greedy/FP-greedy) all 50 trials are bit-identical — those cells are N=1. `results.md:12` claims "independent initial conditions." | **OPEN** — Gate B (B1). |
| D3 | **Candidate trajectories may not reach their vertices.** `min_Δv_dist` folds terminal miss distance into the *objective* at weight 1.0 instead of constraining it; `sum_norm_Δv` returns Σ‖Δv‖² despite its name. | `min_fns.jl:54-55`, `obj_fns.jl:160-171` | If the miss is comparable to R=6.378 km the action space is fiction and both the lifting proposition and the safety claim collapse. Reviewer 5 #3. **Existential — measure in week 1.** | **RESOLVED — NOT A DEFECT.** Measured 2880 on-policy samples: median miss **4.1 cm**, p95 8.1 cm, max 18 cm = 2.8e-5·R. See note [a]. |
| D4 | **Stage cost is not fuel.** Code computes `0.1*norm(u1-u2)` — the norm of the *difference* of control vectors, which rewards the evader for thrusting differently from the pursuer and the pursuer for matching thrust direction. Paper Eq. (9) writes `λ₂(‖u_P‖ − ‖u_E‖)`, which is differential fuel and *is* meaningful. | `matrix_game_solver.jl:179-190` | The implemented objective has no physical interpretation. Reviewer 7 minor #1. | **FIXED** Aug 31. `λ1·sqrt(dist+0.1) + λ2·(‖u_P‖−‖u_E‖)`; λ1/λ2 in `params`. Tests assert equal burns cancel exactly — the property `norm(u1-u2)` violated. |
| D5 | **`Meta_*` pursuer transpose.** P2's meta branches use `players[2].cost * predicted_v_probs`; the FP branches correctly use `cost'`. | `matrix_game_solver.jl:394-405` | Indexes the wrong axis. Invalidates all existing `Meta_*`-as-P2 data. | **FIXED** Aug 31. Transpose added and `argmin`→`argmax` in `Meta_greedy`; `Meta_mixed` weight formula flipped to match `FP_mixed`. See [e]. |
| D6 | **ΔV constraint inert.** `Δv_max = 2.0` km/s default, never threaded from `params`; observed max per-segment ‖Δv‖ ≈ 0.591 km/s. The paper states 0.1 km/s. | `min_fns.jl:47`, `matrix_game_solver.jl:158` | The constraint never binds, so "comparable ΔV budgets" is unsupported. Four separate literals would need editing to change it. | **CONFIRMED, reframed.** Measured max ‖Δv‖ = **0.0221 km/s**, not 0.591. Plumbing fixed (keyword args + `params.Δv_max`); see [b]. |
| D7 | **Silent LP failure.** The `OPTIMAL` check is a non-fatal `println`; the `error` is commented out. OSQP is a first-order ADMM QP solver (default tolerance ~1e-3) being used on a pure LP, with a `z ≥ 1e-4` floor. | `matrix_game_solver.jl:103-128` | A bad solve is currently invisible, and the weights feed `ProbabilityWeights` and will feed exploitability arithmetic. Severity **unmeasured** — see A1. | **QUANTIFIED and FIXED.** 0.83% `ALMOST_OPTIMAL` under OSQP. Solver swapped to Ipopt; see [c]. |
| D8 | Live `@exfiltrate` in an analysis path; `load_games_vec` reads a path layout `save_games_vec` no longer writes. | `Utils.jl:301`, `plotting.jl:827`, `Utils.jl:367-378` | Analysis drops into Infiltrator; loader is dead code. | **FIXED** Aug 31. `@exfiltrate` removed; `load_games_vec` repaired to the `n<k>/` layout. |
| D9 | **Derivative-free inner solver.** The trajectory optimizer is a hand-rolled augmented Lagrangian whose inner solve is `Optim.optimize(fn, dfn, x_0, NelderMead())` — Nelder-Mead on a 30-dimensional problem (N=10 segments × 3), discarding the ForwardDiff gradient `dfn` that is built and passed to it. | `min_fns.jl:109-124`, `aug_L.jl:184-196` | Nelder-Mead stagnates above ~10 dimensions. Likely a direct cause of D3 (poor terminal miss) and of D6 (the ΔV penalty never activating). Reviewer 7 #1 asked for solver details and real-time suitability. | **FIXED** Aug 31. Replaced by BFGS + BackTracking with analytic AD gradients: 12/12 converged, miss 7.6e-12 km, 38% less fuel. Ipopt was tried and does not converge here — see [h]. |
| D10 | **AD path broken — no gradient exists.** `cart2kep` allocates `oe = zeros(6)` (a `Vector{Float64}`) then assigns orbital elements into it. It sits on the objective's call path via `prop_kepler_tof:125` → `prop_kepler_tof_Nseg:174`, so `ForwardDiff.gradient` of the trajectory objective throws `MethodError` on `setindex!` with a `Dual`. Invisible because the only consumer of that gradient was `NelderMead()`, which never evaluates it. | `propagator.jl:179-186` | **Invalidates §9 risk 4.** No gradient-based solver — Ipopt, LBFGS, or the existing `min_bfgs` — could ever have run. The plan's LBFGS fallback shared the same broken dependency, so it was not a hedge. | **FIXED** Aug 31. One line. Gated: bitwise identical on 20,000 random states, and ForwardDiff now matches `central_fdm(5,1)` to 1.6e-9. |
| D11 | **The fuel objective is inert.** In `min_Δv_dist` the terminal-miss term is a distance in km (1–40 at the Lambert guess) while `sum_norm_Δv` returns Σ‖Δv‖² in (km/s)² (~1e-3). At the implicit weight of 1.0 the miss term outweighs fuel by ~3 orders of magnitude. | `min_fns.jl:74-75`, `obj_fns.jl:159-171` | `min_Δv_dist` does not minimize Δv in any meaningful sense — it is a pure targeting solve. "Comparable ΔV budgets" is unsupported, and D4's stage cost operates on Δv values that were never fuel-optimized. **This is the real defect D3 was standing in front of.** | **FIXED** Aug 31. True ΔV with an ε=1e-12 desingularization for the gradient; `w_miss` in `params`, swept. See [g]. |
| D12 | **Augmented Lagrangian double-counts the objective.** `min_aug_L_eq_ineq` sums `aug_L_fn` and `aug_L_ineq_fn`, each of which already includes `obj_fn(x)`, giving `2·obj_fn(x)` + both penalty sets. The correct single-call form (`aug_L_eq_ineq_fn`) exists and is commented out one line below. | `aug_L.jl:241-243` | The objective is weighted 2× against the constraints, so both are systematically under-enforced. Off the game hot path (`min_Δv_dist` → `min_aug_L_ineq`), but reached by `min_Δv` — the "correct pattern" for hard-constraining miss distance. | **OPEN**, deferred. The AL is now bypassed on the hot path (0/12 escalations), so this is reachable only via `min_Δv`. Fix if that path is ever used. |
| D13 | **Opponent identification is structurally uncomputable for FP opponents.** `tracked_strategies` is `["mixed","greedy","random",1..6]` — it contains no `FP_*` or `Meta_*` hypotheses. | `IC.jl:74` | `p1_correct_id_rate` / `p2_correct_id_rate` are `NaN` against every FP opponent; `test/MC_results_late_game.csv` confirms this for all FP columns. `results.md` reports "correct-ID rate is 1.0 for all identifiable opponents" — true only because FP opponents are excluded by the word *identifiable*. The metric was never computed for the opponents the paper is about. | **OPEN** — a modeling gap, not a bug: extending it needs the opponent's `fp_belief`, which a player does not observe. Surface in the writing. |

### Notes

**[a] D3 is dead — the lifting is real.** Measured Aug 31 over 2880 on-policy samples (8 seeds x
30 steps x 2 players x 6 vertices, `random/random`):

| | median | p95 | max |
|---|---|---|---|
| terminal miss | 4.1 cm | 8.1 cm | 18 cm |
| as a fraction of R = 6.378 km | 6.5e-6 | 1.3e-5 | 2.8e-5 |

No late-game degradation (steps 25-30 indistinguishable from 1-6), no player asymmetry. The
decision rule `p95 < R/2` passes with ~40,000x margin. **Section 9 risk 2 does not fire.** The
miss is now recorded per solve in `player_struct.solve_info.traj`, so it is monitored
continuously rather than by an offline probe, and `test/runtests.jl` asserts it.

**[b] D6 reframed — the ΔV cap is irrelevant, not merely un-threaded.** Measured max per-segment
‖Δv‖ = **0.0221 km/s**. The original claim of 0.591 km/s is **not reproduced** - it is 27x
smaller. Critically, **the paper's stated 0.1 km/s would also never bind**; binding requires
~0.02. So the honest fix is not "set it to 0.1 and now the constraint is real" - it is either to
state that the constraint is inactive at this scenario scale, or to choose a cap reflecting an
actual thruster. Plumbing is fixed either way: `dm`/`Δv_max` are now keyword-only, and
`params.Δv_max` / `params.w_miss` thread to the call site.

**[c] D7 measured, and the solver replaced.** 0.83% of OSQP solves returned `ALMOST_OPTIMAL`,
silently, into `ProbabilityWeights`. A solver-independent `nash_certificate` (saddle-point
residual; consults no status flag) was added and the backends benchmarked over 300 random 6x6
games:

| backend | median residual | fraction > 1e-8 | ms |
|---|---|---|---|
| OSQP, defaults (**the CDC-era solver**) | 1.7e-3 | **100%** | 22.8 |
| OSQP, polish + eps 1e-9 | 2.6e-16 | **24%** (bimodal) | 9.4 |
| OSQP, polish + eps 1e-12 | 2.6e-16 | 25% (bimodal) | 21.4 |
| **Ipopt, tol 1e-14** | **2.7e-14** | **0%** | **3.3** |

OSQP-with-polish is exact when polish succeeds and ~5e-3 when it does not; tightening tolerances
does not move that split. **Every CDC-era LP was solved to roughly three decimal places.** Ipopt
is uniformly accurate, fastest, and already a dependency - no Manifest change. This table is the
solver characterization Reviewer 7 #1 asked for. Failure policy is now record-and-fall-back
(never `error()` inside `Threads.@threads`), so the failure *rate* becomes a reportable number.

**[d] A4 sign trap.** The "intended" form commented at `matrix_game_solver.jl:324` is
`0.1*(norm(u1) - norm(u2))`. Under the corrected D1 orientation `stage_cost` is P1's payoff and
P1 *maximizes* it, so with `u1 = u_E`, `u2 = u_P` the paper's Eq. (9) requires
`lambda2*(norm(u2) - norm(u1))` - the commented line is **sign-inverted**. Uncommenting it would
survive A1 and quietly poison Gate B. **Fixed Aug 31**, and pinned by the "stage cost
orientation" testset, which asserts the property the old form actually violated: two equal burns
must cancel exactly, including equal-magnitude burns in orthogonal directions.

**[e] D5 is three defects, not two.** P2's `Meta_greedy` (`:562-566`) is missing the transpose
*and* uses `argmin` where the correct `FP_greedy` P2 branch (`:550-553`) uses `argmax` - since
`players[2].cost = -A`, `argmax` is right. `Meta_mixed` (`:568-573`) is missing the transpose
*and* has its weight formula inverted relative to `FP_mixed`. Because both dimensions are 6, the
axis error produces a plausible-looking index instead of erroring.

**[f] The game is NOT degenerate — but the decision margin is thin, and there is a geometric
confound.** Measured Aug 31 over 240 matrices (4 matchups x 2 seeds x 30 steps), normalized by
`|mean(A)|`:

| metric | median | p90 | max |
|---|---|---|---|
| full matrix range | 46.0% | 59.3% | 88.2% |
| P1 leverage vs a uniform opponent | 15.4% | 21.7% | 25.1% |
| **best-vs-2nd-best row gap** | **1.7%** | 5.1% | 7.4% |

Leverage grows from 8.3% (steps 1-6) to ~16-17% from step 7 on, consistent across all four
matchups. And per-step edges **compound**: pure-vertex strategies (P1 always plays vertex j) vs a
`mixed` pursuer give last-8-step separations spanning **39.5% and 47.2%** of the mean across the
six vertices, on two seeds. So strategy selection matters a great deal over 30 steps.

Two consequences:

1. **The 1.7% decision margin is the number to worry about.** The D1 sign error produced
   NashConv ~0.1 absolute on values of ~2.3 -- about **2.5x the decision gap** -- so it routinely
   flipped which vertex was chosen. OSQP's 1.7e-3 residual is ~15-20% of the gap: not dominant,
   but not ignorable for `greedy`, which is a bare `argmax` over near-ties.

2. **NEW RISK - a geometric confound for the headline.** The vertex ranking is partly
   seed-independent: vertices 3 (`botin`) and 6 (`topout`) are top-2 on both seeds, vertex 1
   (`top`) is bottom-2 on both. 3 and 6 are the antipodal `axis_2`-dominant pair, i.e. in-plane;
   vertex 1 is pure +`axis_3`, out-of-plane. That is physically expected -- an in-plane radial
   offset converts into **along-track drift that grows over an orbit**, while out-of-plane
   displacement merely oscillates. So the six vertices are not equally valuable for evasion, and
   the asymmetry is static geometry, independent of the opponent.

   **If FP-greedy wins partly by converging onto vertex 3 or 6, some of its measured advantage is
   discovering a fixed geometric bias, not modeling an opponent** -- which is exactly the
   headline claim. An AAMAS reviewer will ask this.

   **Required control (add to Gate B, alongside B6):** a *best-fixed-vertex* baseline -- choose
   the single best vertex in hindsight and play it every step. If FP does not clearly beat it,
   the opponent-modeling claim is in trouble. If it does, the margin over that baseline *is* the
   opponent-modeling effect, cleanly separated from geometry. `choose_strategies!` already
   accepts an `Int` strategy, so this costs nothing to implement.

**[g] D11 fixed — fuel is now genuinely part of the objective.** `sum_norm_Δv` returns true ΔV
(with an ε=1e-12 desingularization so the gradient survives a zero-Δv segment, which plain
`norm` would NaN). Measured at the 12 step-1 subproblems, sweeping `w_miss`:

| w_miss | median fuel (km/s) | median miss (km) | p95 miss / R |
|---|---|---|---|
| 1 (default) | 0.04758 | 7.0e-5 | 2.0e-5 |
| 10 | 0.06206 | 5.6e-5 | 1.4e-5 |
| 100 | 0.06249 | 4.5e-5 | 1.1e-5 |
| 1000 | 0.06825 | 2.5e-5 | 6.2e-6 |

A clean monotone tradeoff, and **A3 passes at every setting** by 4+ orders of magnitude.

Two things worth putting in the paper:

- **The optimized trajectory costs *more* fuel than the Lambert initial guess** — 0.0476 vs
  0.0293 km/s, i.e. 162%. That is not a regression: the Lambert guess misses its vertex by a
  median of **20.3 km**, and the extra 63% ΔV is the price of actually arriving. This is a clean
  one-line justification for why the optimizer exists at all, and it retires the D9 worry that
  Nelder-Mead might not be earning its keep on the targeting objective.
- **The balance has flipped.** Before: fuel ~1e-3 (squared norms) against a miss starting at
  1–40 km, so miss outweighed fuel ~1000:1 and ΔV was not optimized. Now fuel ~0.048 against a
  converged miss of 7e-5, so both terms are active — miss dominates early (driving the solve to
  the vertex) and fuel dominates near the solution (trimming waste).

**Open question for the author, not a defect:** `w_miss` could go *below* 1 to buy back fuel,
since the miss has four orders of magnitude of headroom before the A3 threshold. Left at 1.0
because the action semantics want the vertex actually reached; worth a sentence either way.

**[h] A6 resolved — but with quasi-Newton, not Ipopt.** Ipopt was tried first and **does not
converge on this problem** at any iteration budget up to 5000 (the terminal miss even oscillates
non-monotonically: 1e-6 at 200 iterations, 4e-3 at 1000, 2.9e-5 at 5000). Cause: the problem is
stiff — a 1e-5 km/s change in Δv moves the terminal position by metres over 10 segments — and a
JuMP user-defined operator gives Ipopt no structure, so `hessian_approximation="limited-memory"`
cannot build a usable curvature model. Measured on the 12 step-1 subproblems (median):

| solver | converged | wall | terminal miss | ΔV |
|---|---|---|---|---|
| AL + Nelder-Mead (incumbent) | n/a | 0.122 s | 7.0e-05 km | 0.04758 |
| **BFGS + BackTracking** | **12/12** | 0.154 s | **7.6e-12 km** | **0.02926** |
| LBFGS + BackTracking | 12/12 | **0.016 s** | 3.2e-12 km | 0.04169 |
| Ipopt (JuMP `@operator`, constrained) | **0/12** | 0.162 s | 1.0e-06 km | 0.03187 |

**BFGS is now the default**: 38% less fuel and ~10^7 better terminal miss than the incumbent, for
26% more wall time (~36 min per 5x5 sweep). Two details were load-bearing, and neither is in the
original plan:

- **A guarded objective.** A large enough trial Δv makes the Kepler propagation non-finite, and a
  NaN trips an assertion *inside* Optim's line search (`isfinite(phi_c)`). Returning a large
  finite value lets the line search back off.
- **`BackTracking`, not the `HagerZhang` default.** HagerZhang extrapolates into the blow-up
  region; backtracking only shrinks the step. This alone moved convergence from 1/12 to 12/12.

Also: **the augmented Lagrangian is now bypassed unless the ΔV cap actually binds** (verified
post-solve; 0/12 escalations, consistent with note [b]), and `min_aug_L_ineq` gained the
50-iteration cap it never had. `min_Δv_dist_solve` returns `converged`/`iters`/`escalated` per
solve, recorded in `solve_info.traj`, which is what makes a sweep-wide convergence rate
reportable.

**For the paper**, the honest and stronger claim is the comparison itself: *"we evaluated an
interior-point NLP solver (Ipopt) via a user-defined-operator interface; it failed to converge on
this problem class at any iteration budget, while quasi-Newton with analytic AD gradients and a
backtracking line search converged on 100% of solves."* That answers Reviewer 7 #1 better than
either result alone. Caveat to state: `g_converged` is 0/12 — these terminate on step/objective
tolerance, not gradient norm, because the terminal-miss term is an **exact penalty** (a norm, so
non-smooth precisely at the solution). A squared miss is smooth but strictly worse in practice —
its gradient vanishes near zero, so the optimizer stops driving the miss down (measured 4.2e-2
instead of 1.4e-9).

**[i] Smoke sweep (Day 9) — the pipeline works end to end.** 25 matchups x 5 games x 10 steps
= 125 games, **15,000 trajectory solves, 2,500 LP solves, 27.3 min** on 32 threads.

| trajectory optimizer (BFGS + BackTracking) | |
|---|---|
| converged | **99.53%** (14,929 / 15,000) |
| `g_converged` | 0.01% — terminates on step/objective tol, see [h] |
| escalated to augmented Lagrangian | **0.00%** — the ΔV cap never binds, confirming [b] |
| iterations | median 820, p95 1,303, max 2,000 (the cap) |
| terminal miss | median 7.3e-12 km, p95 6.1e-10, **max 8.2e-05** (R = 6.378) |
| total ΔV per solve | median 0.0102 km/s |

| matrix-game LP (Ipopt, tol 1e-14) | |
|---|---|
| certificate violation | median 3.7e-15, **max 1.1e-12** (threshold 1e-8) |
| above threshold | **0 / 2,500** |
| fell back to pure maximin/minimax | **0 / 2,500** |
| maximin − minimax gap | max 7.5e-13 |
| wall per LP | median 0.0057 s |

**`mixed`/`mixed` NashConv of the strategies actually played: median 1.6e-14, max 3.3e-12.**
The played mixed strategies *are* the stage-game equilibrium — D1 is fixed end to end, not just
in unit tests. The 99.53% convergence rate is the number Reviewer 7 #1 asked for.

**COST WARNING for Gate B.** In-pipeline median wall per trajectory solve is **0.505 s**, not the
0.154 s measured on step-1 subproblems — later steps are harder (median 820 BFGS iterations, and
some hit the 2,000 cap). Projecting a full 5x5 at 50 games x 30 steps on 32 threads:

    per game    30 steps x 12 solves x 0.505 s  = 182 s
    per matchup 50 games over 32 threads (2 waves) = ~364 s
    full sweep  x 25 matchups                    = ~2.5 hours

That is ~5x the earlier estimate and ~5x the Nelder-Mead incumbent. Levers, in order of
preference: **LBFGS** instead of BFGS (measured 0.016 s vs 0.154 s on step-1 subproblems, ~10x
faster, at ~40% more ΔV — the accuracy is unaffected, both reach ~1e-12); lower `maxiter` from
2,000; or accept 2.5 h, which is still an overnight-free single sitting. Decide before B3.

**Not a defect, but state it in the paper.** `axis_123` is *not* an orthogonal frame:
`a1 . a2 = -9.9e-3 ~ -e`. The hexagon lies in the `a2`-`a3` (radial/normal) plane, **not** the
plane normal to velocity; those coincide only for a circular orbit and differ by ~0.57 deg at
e=0.01. Reviewer 5 #1 asked exactly how vertex positions are obtained, so say which plane.

---

## 4. Workstreams

### Gate A — correctness + existential diagnostics (Aug 31 – Sep 6)

Runs on **existing** `.jld2` files. No re-run needed yet.

**A1. Fix the LP (D1), instrument the solver (D7).**
- `solve_mixed_nash(A)` → `sms(-A)` for the evader (maximizer of A), `sms(A')` for the pursuer (minimizer of A). Return the value, currently discarded at `:58`.
- **Instrument before swapping solvers.** Log LP residual, `sum(w)`, `min(w)`, and termination status across a full run. Restore the `error(...)` at `:125` regardless. Swap `OSQP.Optimizer` → `HiGHS.Optimizer` (`:107`, one line, add to `Project.toml`) **only if the logs justify it** — but report the solver characterization either way, since Reviewer 7 #1 asked for exactly this.
- Re-examine the `z[i] >= 1e-4` floor (`:118-120`): it distorts the equilibrium and the paper describes it as a "probability floor," which it is not. Drop it or state it honestly as ε-support regularization.
- **Regression test — the one that would have caught D1.** Von Neumann: maximin = minimax. Assert `-sms(-A).v ≈ sms(A').v ≈ p'Aq` on random matrices and on saved `A_k`. `test/runtests.jl` is currently 6 lines and tests nothing.

**A2. Exploitability machinery.** New `src/game_theory/exploitability.jl`. For evader payoff `A`, evader mix `p`, pursuer mix `q`: evader BR gain `max_i(Aq)_i − p'Aq`; pursuer BR gain `p'Aq − min_j(p'A)_j`; `NashConv` = sum; exploitability = NashConv/2.
- Everything needed is already persisted per step: `game.p1_state[k].cost` = `A_k`, `.weights` = `p_k`, `game.p2_state[k].weights` = `q_k`, `.fp_belief` = empirical opponent frequency. **Renormalize `weights` first** — they are raw LP output.
- Run against existing `test/results/n30/` to quantify what D1 cost. Internal number; will not appear in the paper.

**A3. Miss-distance diagnostic (D3) — the existential check. Run this twice: before A6, and again after.** Recompute vertices post-hoc via `polygon_vertices(game.rv_ref_E[k][end,:], params)`, compare against `p.X[j][end,1:3]` for every step/player/vertex. Report median and 95th percentile as a fraction of `R_polygon` = 6.378 km. *Verify the `rv_ref_E[k]` ↔ `p1_state[k]` index pairing first — `prop_game_step` pushes the reference before the player states.*
- **Pass 1 (existing data):** establishes the status quo under Nelder-Mead.
- **Pass 2 (after A6):** the Ipopt switch is expected to reduce terminal miss substantially, and may resolve D3 outright.
- **Decision point.** Median miss ≪ R → proceed. Comparable to or greater than R *after* A6 → the hexagon does not describe the action space, and the terminal miss must become a hard constraint rather than an objective term.

**A4. Fix the stage cost (D4).** Change to the paper's differential-fuel form `λ₂(‖u_P‖ − ‖u_E‖)`. **Keep fuel in the objective** — it is physically essential in this domain. Report ΔV separately as a metric, and handle Reviewer 7's "partly overlapping metrics" point in the writing, not the model.

**A5. Housekeeping + runtime.** Remove live `@exfiltrate`s (D8); repair or delete `load_games_vec`; add `@elapsed` around the 12 `min_Δv_dist` calls per step and around the LP. There is currently **zero** timing instrumentation anywhere in the repo, and Reviewer 7 #1 asked for runtime characterization. Instrument *before* A6 so the switch can be quantified.

**A6. Replace Nelder-Mead with Ipopt (D9). Must land before the Gate B re-run.**

Ipopt is already a dependency (v1.14.0 in `Manifest.toml`) and `src/Opt/min_ipopt.jl` exists but has zero call sites in the game pipeline. The intended restructure is to **drop the hand-rolled augmented Lagrangian entirely** rather than keep it with a different inner solver — Ipopt is an interior-point NLP solver that handles inequality constraints natively:

- **Variables:** `x ∈ R^{3N}`, N=10 segments (30 decision variables).
- **Objective:** `sum_norm_Δv(x,N) + miss_distance_prop_kepler_Nseg(rv_0, x, N, rv_f, tof_N, mu)`.
- **Constraints:** `‖Δv_i‖₂ ≤ Δv_max`, i = 1..N — ten nonlinear inequality constraints, enforced properly instead of via a penalty that currently never activates.

This directly answers Reviewer 7 #1, which asked whether the MPC problem is convex and reliably solvable. **Report honestly: it is a nonconvex NLP solved to local optimality** (Kepler propagation makes the dynamics nonlinear). With Ipopt you can report per-solve convergence status, iteration counts, wall time, and — importantly — the **fraction of the ~18,000 solves per matchup that converge**, plus the fallback behavior on failure.

*Derivative plumbing is the main implementation risk.* The objective propagates through `prop_kepler_tof_Nseg`, which contains an iterative Kepler solve. It is already ForwardDiff-compatible — `min_fns.jl:117` builds a ForwardDiff gradient of exactly this objective today — so AD through it works; the question is only whether it survives JuMP's nonlinear interface. Register the objective as a user-defined operator with a ForwardDiff gradient, or go through `ADNLPModels`. **Fallback if the JuMP/Ipopt AD path fights back:** keep the augmented Lagrangian but swap the inner `NelderMead()` for `LBFGS()` in `min_fns.jl:112`, which still discards Nelder-Mead and still uses the gradient, at a fraction of the effort.

**Consequences to propagate:** every trajectory changes, so all existing results are invalidated — acceptable, since D1 already invalidates them. **Re-measure per-matchup wall clock after this lands**, because the Gate B and Gate C sweep estimates in §1 and §8 all derive from the 6 min/matchup Nelder-Mead baseline. Expect it to drop. Also revisit D6: with the constraint enforced natively, a physically meaningful `Δv_max` becomes settable and can actually bind.

### Gate B — the honest re-run (Sep 7 – Sep 13). **Narrative locks at the end of this gate.**

**B1. Real Monte Carlo (D2).** Add IC dispersion as keyword args to `init_game`, driven by the `local_rng` already passed at `play_games.jl:104`. Do **not** simply uncomment `rand_IC` — it drops the `a_P = 1.005a` drift and hands the pursuer the evader's velocity, a different scenario. Disperse around the nominal instead: pursuer semi-major-axis ratio over a range, plus phase/RAAN jitter on both craft. Raise trials to 100–200. (Peters et al. use 50–100 randomized initial states with SEM ribbons — that is the bar.)

**B2. Slim persistence.** Files are 59 MB each. At 200 trials × 49 matchups the full candidate-trajectory dump is untenable. Add a slim save retaining `A_k`, both weight vectors, `chosen`, beliefs, executed states, and terminal misses — everything the analysis actually reads.

**B3. Re-run the 5×5 with correct signs**, then compute exploitability, NashConv, and capture rate (post-hoc: fraction of games whose minimum inter-player distance falls below a stated capture radius — there is currently **no capture or termination condition anywhere in the codebase**, and a pursuit-evasion reviewer will ask what capture means).

**B4. Confidence intervals and significance.** `MC_stats` computes std but `print_MC_stats` prints means only, and `analyze_late_game.jl` reduces to `mean` before writing the CSV, so all dispersion is discarded. Add CI columns and a paired test on headline comparisons. **No number enters a table without an interval.**

**B5. EGTA meta-game.** Treat the strategy table as a normal-form meta-game; report its Nash equilibrium (optionally α-Rank). Pure post-processing, and the principled fix for Reviewer 8's objection that "dominant strategy" is misused — afterwards you can say precisely whether FP-greedy is dominant or merely in the support.

**B6b. Best-fixed-vertex baseline (added Aug 31 — see note [f]).** Play a single hindsight-best vertex every step. Controls for the static geometric asymmetry between hexagon vertices, without which the opponent-modeling headline is confounded. Uses the existing `Int` strategy branch.

**B6. Oracle best-response baseline.** A player told the opponent's true strategy that best-responds exactly. ~10 lines, and it establishes the ceiling: "FP recovers X% of the oracle's advantage within N steps." This is the cheap stand-in for a MARL baseline (see §5).

### Gate C — the contribution (Sep 14 – Sep 20)

**C1. Recency-weighted FP.** `f_s = γ·f_{s-1} + g_{s-1}` in `update_beliefs!` (`matrix_game_solver.jl:444-455`), γ in `params`, γ=1 recovering current behavior; effective memory ≈ 1/(1−γ). Sweep γ. The argument: uniform averaging pools observations across game steps with *different* payoff matrices — at step 30 the players are three orbits from step 1.

**C2. A no-regret learner.** Regret matching (`p_j ∝ [cumulative regret_j]⁺`) and/or Hedge (`p_j ∝ exp(η · cumulative payoff_j)`) as new branches in `choose_strategies!`. Rationale: it is the baseline AAMAS expects; in zero-sum games two no-regret players' time-averaged play converges to Nash at O(1/√T) where FP's worst case is far worse; and **dynamic-regret** bounds exist for time-varying games, which is what makes C3's tracking claim a theorem rather than an observation. Requires a persistent learner-state field on `player_struct`, which uses positional construction — `IC.jl:76` and `play_games.jl:27-28` must be updated together.

**C3. FP convergence verification.** Two parts; the distinction is the point.
- *Frozen-matrix sanity check.* Freeze `A_k` at a representative step, run FP for T iterations, show ‖q̂_T − w*‖ → 0 and exploitability → 0. Validates the FP implementation and the LP jointly. **Robinson (1951) guarantees convergence only for a fixed payoff matrix** — this is the only setting where the classical theorem applies.
- *Time-varying tracking.* Measure ‖A_{s+1} − A_s‖ to quantify how slowly the game varies, then plot per-step exploitability of empirical play against it. The honest claim is a **tracking** result: empirical play follows the moving stage-game equilibrium with a gap bounded by the variation budget. **Do not claim Robinson applies here.**
- State explicitly the assumption that makes cross-step pooling meaningful at all: vertex *j* denotes a consistent direction in the local ê₂–ê₃ frame across steps even though it is a different physical point each step. Reviewers will look for this, and it is the honest weak point of the method.

**C4. Vertex ablation, 6 / 12 / 24.** Four hardcoded sites: `polygon_vertices` (`Utils.jl:100-135`) returns a fixed-field NamedTuple `(top, topin, botin, bot, botout, topout)`, so arity is baked into the return type; `zeros(6,6)` (`matrix_game_solver.jl:207-208`) even though `n_vertices` is computed correctly at `:210`; `fp_belief = ones(6)` (`IC.jl:76`); fixed-vertex hypotheses `1..6` in `tracked_strategies` (`IC.jl:74`). Everything else derives arity from `size`/`length`. Cost scales with vertex count: ~6 h at 12, ~12 h at 24 — overnight jobs. Purpose is a convergence curve showing the lifted value stabilizes as the action set refines, so the hexagon is not arbitrary; 6/12/18 is an acceptable fallback.

**C5. Lifting proposition.** Restricting a maximizer's action set can only lower the value (a minimizer's can only raise it) — a one-line rigorous bound. A Lipschitz stage cost over a δ-net of the reachable set gives |V_M − V_cont| ≤ L·h_M with h_M the covering radius. C4 supplies the empirical curve. Together, the direct answer to Reviewer 5 #4.

### Gate D — writing support (Sep 21 – Oct 8)

Deliver to Overleaf as markdown prose, CSV tables with intervals, and PDF figures. Optional if time allows: promote `Meta_greedy`/`Meta_mixed` from "future work" to a full type-based-reasoning comparison (requires D5 first) — squarely AAMAS's territory via Albrecht & Stone, but it costs page space in 8 pages and expands the sweep to 7×7 ≈ 6 h.

**Double-blind logistics:** the repo is public and identifiable, and `RL_BENCHMARK_HANDOFF.md` names the author. Supplementary material (≤25 MB zip) must be anonymized — author names and git history stripped.

---

## 5. Decisions already made — do not relitigate

- **MARL baseline: no, not for this deadline.** Each game step costs 12 nonlinear trajectory optimizations; RL needs 10⁵–10⁶ environment steps, which would require a surrogate environment that is then not this problem. It is also an unfair comparison as usually run — an offline-trained state-conditioned policy against an online adapter given no training. Use the oracle best-response (B6) instead. MARL is the direction `RL_BENCHMARK_HANDOFF.md` already points, as separate work.
- **Fuel stays in the objective.** Fix its functional form (D4); do not remove it.
- **Nelder-Mead is out; Ipopt is in (A6).** Decided by the author, not gated on measurement. Fallback is Optim + LBFGS, never Nelder-Mead.
- **Do not claim time-variation as novelty.** See §0.
- **Venue is AAMAS**, not an aerospace conference. Decided by the author.

## 6. Drop order if the calendar slips

Cut from the bottom: **`Meta_*` promotion** → **24-vertex ablation** (keep 6/12) → **no-regret learner** (keep discounted FP, a one-line change) → **EGTA meta-game**.

**Never cut:** the sign fix (A1), the miss-distance diagnostic (A3), randomized ICs with confidence intervals (B1/B4), exploitability (A2). Those four are what make the headline provable rather than asserted, and the Monte Carlo defect in particular is the kind of thing that sinks a paper when a reviewer finds it.

## 7. Verification

- `test/runtests.jl` gains real tests: maximin = minimax on random and saved matrices; weights sum to 1 and are non-negative; solver termination is `OPTIMAL`; zero-sum property `p2.cost == -p1.cost`; exploitability of the LP solution ≈ 0 by construction (the strongest single check on A1).
- Frozen-matrix FP convergence (C3) doubles as an integration test: FP's empirical frequency must converge to the LP's `w*` on a fixed matrix. Under current code these disagree, so it fails before the fix and passes after.
- Re-run one previously-computed matchup with the fix disabled and confirm the old numbers reproduce, establishing the pipeline is otherwise unchanged.
- Spot-check terminal miss distances against `R_polygon` before and after any change to the trajectory optimizer.
- Log `-sms(-A).v` against `sms(A').v` at every step of a full run; flag any divergence.

## 8. Schedule

| Window | Work | Gate |
|---|---|---|
| Aug 31 – Sep 6 | Task zero (`CLAUDE.md`); A1–A5 | **Is the lifting real?** |
| Sep 7 – Sep 13 | B1–B6 | **Narrative locks** |
| Sep 14 – Sep 20 | C1–C5 | OpenReview accounts by ~Sep 17 |
| Sep 21 – Sep 30 | Writing support; figures; abstract drafted | |
| **Oct 1** | **Abstract deadline** (100–300 words; title and authors lock) | |
| Oct 1 – Oct 8 | Polish, anonymized supplementary, buffer | **Oct 8 paper deadline** |

## 9. Risks

1. **The headline may shrink.** With correct signs, `mixed` becomes genuinely unexploitable and FP's gain against it should approach zero — that is what minimax *means*. The exploitation-vs-safety framing is robust to this because it predicts it, but the abstract cannot be written before Gate B closes.
2. ~~**D3 could be fatal.**~~ **RETIRED Aug 31.** Measured: median terminal miss 4.1 cm against R = 6.378 km, no late-game degradation. The action space is well defined and the lifting is faithful. See note [a].
3. **Novelty pressure from [3].** Peters et al. covers lifting, mixed strategies, receding horizon, time-varying payoffs, and a randomized-IC tournament with SEM. The delta is online opponent modeling and the exploitation/safety tradeoff. Everything in the paper should be pointed at that delta; anything that reads as "Peters et al. applied to orbits" should be cut.
4. **The Ipopt migration (A6) is the largest single code risk.** It touches the innermost, most-called function in the system and perturbs every trajectory in every result. **One stated mitigation was false and has been corrected:** the objective was *not* ForwardDiff-differentiable (D10), and the LBFGS fallback shared that same broken dependency, so it was never a hedge. Both are fixed as of Aug 31 — AD now matches finite differences to 1.6e-9 — so the mitigation is real now rather than assumed. The remaining mitigations stand: it lands in Gate A before the expensive sweeps, and LBFGS is a small change if JuMP's nonlinear interface resists. If A6 slips past Sep 6, take the fallback rather than pushing the Gate B re-run.
5. **Ipopt may change the story, not just the numbers.** Better convergence should lower terminal miss (helping D3) and make the ΔV constraint bindable (fixing D6) — but it also means the corrected `mixed`/`greedy` strategies are being computed on a *different, better* action set than the CDC results used. Treat every CDC number as void; do not attempt reconciliation.
