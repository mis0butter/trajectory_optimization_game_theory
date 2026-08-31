# AAMAS 2027 submission — handoff plan

**This document is written to be handed to an agent on a different machine.** It assumes no prior context about this repository.

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
- **Disk:** `test/results/` currently holds ~12 GB. Gate A diagnostics need only `test/results/n30/` (~1.4 GB, 25 files at ~59 MB) copied to the new machine. Budget 20–40 GB for new runs unless B2 (slim persistence) lands first.
- **Measured baseline cost:** the 5×5 sweep at 50 games × 30 steps took **2 h 51 min** wall clock on an 8-core machine, ~6 min per matchup, derived from output file mtimes.

---

## 2. Task zero — `CLAUDE.md`

Before any code changes, write `CLAUDE.md` at the repo root (README.md is currently empty and open in the user's editor; write `CLAUDE.md` and leave README alone unless asked). It must capture:

- **Purpose and the game:** evader/pursuer, hexagonal lifting, per-step matrix game, receding horizon.
- **Code organization:** `src/game_theory/` (matrix game solver, play loop), `src/Opt/` (trajectory optimization: augmented Lagrangian, objective/constraint functions), `src/Dyn/` (Kepler propagation), `src/Lambert/` (initial guess), `src/Utils/` (ICs, structs, plotting, MC statistics), `src/old/` (dead). Entry points: `test/run_all_MC_table.jl` produces the sweep, `test/analyze_late_game.jl` post-processes into `test/MC_results_late_game.csv` and figures.
- **Execution path for one game step**, which is the thing a new agent most needs: `prop_game_step` → `compute_states_nash` → `compute_players_XU` (12 trajectory optimizations) → `compute_cost_matrices` → `compute_mixing_weights!` (LP) → `choose_strategies!` → `update_beliefs!` → `update_chosen_trajectories!`.
- **CDC rejection and reviewer feedback**, summarized, pointing at `Game_Theory_Space_CDC_2026/reviews/`.
- **The defect table from §3 below**, verbatim — this is the highest-value content for a new agent.
- **Paper-source situation:** the LaTeX in `Game_Theory_Space_CDC_2026/` (`root.tex`, `sections/`, `refs.bib`) is the *old pre-FP draft*, not the CDC submission. The submitted version exists in this repo only as `Game_Theory_Space_ICRA_2026 (2).pdf`; its source lives in Overleaf. `test/results.md` is a third, intermediate write-up.
- **Gotchas:** threads default to 1; relative paths; `Manifest.toml` has no LP solver but OSQP; there is no real test suite (`test/runtests.jl` is 6 lines); `test/` is scripts, not tests.

---

## 3. Verified defects

| # | Defect | Location | Consequence |
|---|---|---|---|
| D1 | **LP sign inversion.** `solve_mixed_security_strategy` returns the row player's *minimizing* security strategy. `solve_mixed_nash(A)` therefore gives the evader a distance-minimizing mixed strategy and the pursuer a distance-maximizing one — both reversed. Verified numerically on a matrix with a known saddle point: the evader is assigned ~1.0 weight on its worst row. | `matrix_game_solver.jl:55-59` | Every `mixed` and `greedy` number in the CDC tables is the equilibrium of the reversed game. Reviewer 7 #2, Reviewer 8 #1. |
| D2 | **No Monte Carlo.** `init_game` hardcodes both orbits; `rand_IC` exists but its only call site is commented out. Only the vertex-sampling RNG varies between trials. | `IC.jl:19-26`, `IC.jl:24-25` | For deterministic matchups (greedy/greedy, greedy/FP-greedy, FP-greedy/FP-greedy) all 50 trials are bit-identical — those cells are N=1. `results.md:12` claims "independent initial conditions." |
| D3 | **Candidate trajectories may not reach their vertices.** `min_Δv_dist` folds terminal miss distance into the *objective* at weight 1.0 instead of constraining it; `sum_norm_Δv` returns Σ‖Δv‖² despite its name. | `min_fns.jl:54-55`, `obj_fns.jl:160-171` | If the miss is comparable to R=6.378 km the action space is fiction and both the lifting proposition and the safety claim collapse. Reviewer 5 #3. **Existential — measure in week 1.** |
| D4 | **Stage cost is not fuel.** Code computes `0.1*norm(u1-u2)` — the norm of the *difference* of control vectors, which rewards the evader for thrusting differently from the pursuer and the pursuer for matching thrust direction. Paper Eq. (9) writes `λ₂(‖u_P‖ − ‖u_E‖)`, which is differential fuel and *is* meaningful. | `matrix_game_solver.jl:179-190` | The implemented objective has no physical interpretation. Reviewer 7 minor #1. |
| D5 | **`Meta_*` pursuer transpose.** P2's meta branches use `players[2].cost * predicted_v_probs`; the FP branches correctly use `cost'`. | `matrix_game_solver.jl:394-405` | Indexes the wrong axis. Invalidates all existing `Meta_*`-as-P2 data. |
| D6 | **ΔV constraint inert.** `Δv_max = 2.0` km/s default, never threaded from `params`; observed max per-segment ‖Δv‖ ≈ 0.591 km/s. The paper states 0.1 km/s. | `min_fns.jl:47`, `matrix_game_solver.jl:158` | The constraint never binds, so "comparable ΔV budgets" is unsupported. Four separate literals would need editing to change it. |
| D7 | **Silent LP failure.** The `OPTIMAL` check is a non-fatal `println`; the `error` is commented out. OSQP is a first-order ADMM QP solver (default tolerance ~1e-3) being used on a pure LP, with a `z ≥ 1e-4` floor. | `matrix_game_solver.jl:103-128` | A bad solve is currently invisible, and the weights feed `ProbabilityWeights` and will feed exploitability arithmetic. Severity **unmeasured** — see A1. |
| D8 | Live `@exfiltrate` in an analysis path; `load_games_vec` reads a path layout `save_games_vec` no longer writes. | `Utils.jl:301`, `plotting.jl:827`, `Utils.jl:367-378` | Analysis drops into Infiltrator; loader is dead code. |
| D9 | **Derivative-free inner solver.** The trajectory optimizer is a hand-rolled augmented Lagrangian whose inner solve is `Optim.optimize(fn, dfn, x_0, NelderMead())` — Nelder-Mead on a 30-dimensional problem (N=10 segments × 3), discarding the ForwardDiff gradient `dfn` that is built and passed to it. | `min_fns.jl:109-124`, `aug_L.jl:184-196` | Nelder-Mead stagnates above ~10 dimensions. Likely a direct cause of D3 (poor terminal miss) and of D6 (the ΔV penalty never activating). Reviewer 7 #1 asked for solver details and real-time suitability. |

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
2. **D3 could be fatal.** If candidates do not reach their vertices, the action space is ill-defined and Gate A expands into fixing the trajectory optimizer. Hence the week-1 placement.
3. **Novelty pressure from [3].** Peters et al. covers lifting, mixed strategies, receding horizon, time-varying payoffs, and a randomized-IC tournament with SEM. The delta is online opponent modeling and the exploitation/safety tradeoff. Everything in the paper should be pointed at that delta; anything that reads as "Peters et al. applied to orbits" should be cut.
4. **The Ipopt migration (A6) is the largest single code risk.** It touches the innermost, most-called function in the system and perturbs every trajectory in every result. Mitigations: it lands in Gate A, before the expensive sweeps, so nothing is re-run twice; the objective is already ForwardDiff-differentiable, so the AD path is de-risked; and the LBFGS fallback is a two-line change if JuMP's nonlinear interface resists. If A6 slips past Sep 6, take the fallback rather than pushing the Gate B re-run.
5. **Ipopt may change the story, not just the numbers.** Better convergence should lower terminal miss (helping D3) and make the ΔV constraint bindable (fixing D6) — but it also means the corrected `mixed`/`greedy` strategies are being computed on a *different, better* action set than the CDC results used. Treat every CDC number as void; do not attempt reconciliation.
