# CLAUDE.md

Guidance for Claude Code when working in this repository.

## What this is

`trajectory_optimization_game_theory` — a Julia package implementing a **two-spacecraft orbital
pursuit-evasion game** in low Earth orbit.

An **evader (P1)** follows a reference orbit; a **pursuer (P2)** tries to close distance. Both
start near-circular at 620 km altitude (`a = 6998 km`, `e = 0.01`, `i = 20°`, period ~97 min);
the pursuer is offset at `a_P = 1.005·a`.

At each **game step**:

1. A hexagon of 6 vertices, radius `R_polygon = 6.378 km`, is placed around the evader's
   reference orbit position one horizon ahead.
2. **Both players** solve 6 trajectory optimizations each (12 total) — impulsive Δv sequences
   targeting each vertex.
3. The 6×6 **zero-sum matrix game** `A[i,j]` = mean stage cost with P1 on vertex `i`, P2 on `j`.
4. The matrix game is solved by LP for a mixed Nash equilibrium.
5. Each player picks a vertex per its strategy (`mixed` / `greedy` / `random` / `FP_greedy` /
   `FP_mixed` / `Meta_greedy` / `Meta_mixed`, or an `Int` for a fixed pure vertex).
6. **Receding horizon:** 10 segments are planned, 5 are executed.

This is the "lifting" of a continuous trajectory game onto a finite action set.

**Target venue: AAMAS 2027** (paper deadline Oct 8 2026). Previously submitted to CDC 2026 and
rejected. **`README_plan.md` is the authoritative plan and status document** — read it first.

## Running things

```bash
julia -t 32 --project=.        # ALWAYS pass -t; Threads.nthreads() defaults to 1
```

- **All result paths are relative** — run from the repo root or paths break.
- `run_MC_games_parallel` uses `Threads.@threads`; without `-t` every sweep is serial.
- `test/runtests.jl` — the real test suite (~3080 assertions). Run it after any solver change.
- `test/run_all_MC_table.jl` — the 5×5 matchup sweep, writes `.jld2` under `test/results/`.
- `test/analyze_late_game.jl` — post-processes into `test/MC_results_late_game.csv` + figures.
- Everything else in `test/` is a REPL-style script, not a test.

**Measured cost** (32-core machine): ~56 ms per trajectory solve in steady state, ~0.67 s per
game step, ~46 s per 30-step game single-threaded, **~45 min for a full 5×5 sweep**.

## Code organization

| Path | Role |
|---|---|
| `src/game_theory/matrix_game_solver.jl` | The core. LP Nash solver, `stage_cost`, cost matrices, strategy selection, belief updates |
| `src/game_theory/exploitability.jl` | Exploitability / NashConv / belief exploitability |
| `src/game_theory/play_games.jl` | Game loop: `prop_game_step`, `run_game`, `run_MC_games_parallel` |
| `src/Opt/min_fns.jl` | `min_Δv_dist` (the one the game uses), `min_optim`, `min_bfgs` |
| `src/Opt/aug_L.jl` | Hand-rolled augmented Lagrangian |
| `src/Opt/obj_fns.jl` | `miss_distance_prop_kepler_Nseg`, `sum_norm_Δv`, `constrain_Δv` |
| `src/Opt/min_ipopt.jl` | 27-line unconstrained stub, **not wired into the pipeline** |
| `src/Dyn/kepler.jl` | `prop_kepler_tof_Nseg` — the workhorse propagator |
| `src/Dyn/propagator.jl` | `kep2cart` / `cart2kep`, ODE propagation, `prop_chosen_rv` |
| `src/Dyn/prop_delta_v.jl` | `find_ref_orbit`, `prop_rv_ref` |
| `src/Lambert/` | Battin Lambert solver — supplies the optimizer's initial guess |
| `src/Utils/IC.jl` | `init_game` — **all initial conditions and `params` live here** |
| `src/Utils/structs.jl` | `player_struct`, `game_struct` |
| `src/Utils/Utils.jl` | `polygon_vertices`, save/load, MC statistics |
| `src/Utils/plotting.jl` | GLMakie plotting |
| `src/old/` | **Dead** — not included by the module |

## One game step — the execution path

This is the thing a new agent most needs.

```
prop_game_step                      play_games.jl:3
  find_ref_orbit / prop_rv_ref      prop_delta_v.jl — advance the reference orbit
  prop_chosen_rv                    propagator.jl  — execute last step's chosen trajectory
  compute_states_nash               matrix_game_solver.jl (orchestrator)
    compute_players_XU              12 trajectory optimizations (2 players × 6 vertices)
    compute_cost_matrices           6×6 zero-sum A;  players[2].cost == -players[1].cost
    compute_mixing_weights!         the LP → mixed Nash
    choose_strategies!              per-strategy vertex selection
    update_beliefs!                 FP counts + Bayesian meta-strategy posterior
    update_chosen_trajectories!     truncate horizon (10) to executed (5)
```

Fresh `player_struct`s are built each step; only `fp_belief`, `strategy_belief`,
`tracked_strategies`, and `learner_state` carry forward. `game.p1_state` / `p2_state` are
append-only, so every step's full candidate set is retained — this is why `.jld2` files reach
~59 MB.

## Conventions that will bite you

**Sign convention — the single most important thing in this repo.**

`stage_cost` **increases with separation**. `players[1].cost = +stage_cost`,
`players[2].cost = -stage_cost`. P1 is the evader and wants separation, therefore:

> **P1 (evader) MAXIMIZES `A`. P2 (pursuer) MINIMIZES `A`.**

`solve_mixed_security_strategy(M)` returns the security strategy of the row player who
**MINIMIZES** `M`. So `solve_mixed_nash` uses `sms(-A)` for P1 and `sms(A')` for P2. Getting
this backwards is defect D1, which invalidated every `mixed`/`greedy` number in the CDC
submission. `test/runtests.jl` pins it with a test asserting the *old* pairing is exploitable.

**Other conventions:**

- `player_struct` has an **all-keyword constructor** — use `player_struct(; fp_belief=...)`, not
  the 15-argument positional form.
- The decision vector is a **30×1 `Matrix`** (`reshape(Δv_vec, N*3, 1)`), not a `Vector`, and it
  is **column-major**: segment `i`'s components are at flat indices `i, N+i, 2N+i` — *not*
  `3(i-1)+1 : 3i`.
- `axis_123` is **not an orthogonal frame**: `a1·a2 ≈ -e`. The hexagon lies in the `a2`–`a3`
  (radial/normal) plane, not the plane normal to velocity.
- Vertex order from `polygon_vertices` is `(top, topin, botin, bot, botout, topout)` = 1..6, and
  it returns a **fixed-arity NamedTuple** — changing vertex count means rewriting it (see C4).
- `params` is a plain `NamedTuple` built in `init_game`, stored as a 1-element vector, hence the
  `game.params[1]` idiom.

## Known defects

**`README_plan.md` §3 is the source of truth** — it has locations, consequences, status, and
measured numbers. Compact index:

| # | Summary | Status |
|---|---|---|
| D1 | LP sign inversion — evader given a minimizing strategy | **FIXED** |
| D2 | No Monte Carlo — `init_game` hardcodes both orbits, so deterministic matchups are N=1 | OPEN (Gate B) |
| D3 | "Candidates may not reach their vertices" | **NOT A DEFECT** — median miss 4.1 cm vs R = 6.378 km |
| D4 | Stage cost uses `norm(u1-u2)`, not differential fuel | OPEN (Day 4) |
| D5 | `Meta_*` P2 branches: missing transpose **and** inverted sign | OPEN (Day 8) |
| D6 | ΔV cap inert — measured max ‖Δv‖ 0.0221 km/s vs a 2.0 cap; even 0.1 would not bind | Plumbing FIXED |
| D7 | Silent LP failure; OSQP on a pure LP | **FIXED** — Ipopt, `tol=1e-14` |
| D8 | Live `@exfiltrate`; broken `load_games_vec` | **FIXED** |
| D9 | Nelder-Mead on a 30-D problem, discarding the gradient | OPEN (Day 7) |
| D10 | `cart2kep` allocated `zeros(6)` → **all ForwardDiff broken** | **FIXED** |
| D11 | Fuel objective inert — miss term outweighs fuel ~1000× | OPEN (Day 5) |
| D12 | `min_aug_L_eq_ineq` double-counts the objective | OPEN, deferred |
| D13 | `tracked_strategies` has no FP/Meta hypotheses → ID rate structurally `NaN` | OPEN (modeling gap) |

## Gotchas

- **Threads default to 1.** Always `-t`.
- **Paths are relative.** Run from the repo root.
- **`test/results/` does not exist here** — gitignored, never copied from the original machine.
  Any analysis of "existing results" must regenerate them.
- **The trajectory optimizer is still Nelder-Mead.** Ipopt is currently used *only* for the
  matrix-game LP. Do not confuse the two roles.
- Julia is **1.12.7**; `Manifest.toml` was resolved under **1.11.1**. Any `Pkg.add` re-resolves
  it — treat dependency changes as deliberate, tagged events.
- The repo is on branch **`aamas`**, not `main`. There is no `README.md`.
- `src/old/` and `test/old/` are dead code. `min_golden_ratio` is dead *and* broken (calls `f`,
  parameter is `fn`).
- The Julia LSP reports spurious `UnusedBinding` warnings on NamedTuple fields. Ignore them.

## Paper sources — three divergent write-ups

1. `Game_Theory_Space_CDC_2026/root.tex` + `sections/` — the **old pre-fictitious-play draft**,
   *not* the CDC submission.
2. `Game_Theory_Space_CDC_2026/Game_Theory_Space_ICRA_2026 (2).pdf` — the only in-repo copy of
   the **actual CDC submission**; its source lives in Overleaf.
3. `test/results.md` — a third, intermediate write-up. Several of its claims are now known
   false (notably "independent initial conditions", contradicted by D2).

Reviews are in `Game_Theory_Space_CDC_2026/reviews/` (3 reviewers + AE). Their common theme:
the game-theoretic formulation is under-specified and internally inconsistent. Two of the
complaints turned out to be genuine code bugs (D1, and the solver-detail request behind D7/D9).

**Division of labor: the author owns the paper in Overleaf.** Deliver prose as markdown and
figures as PDF. **Do not create LaTeX in this repo.**

## Working style for this repo

- **Value-preservation discipline.** When fixing infrastructure, verify numbers are unchanged
  (bitwise where possible) so that when results *do* move, the cause is unambiguous.
- **Measure before choosing.** The OSQP→Ipopt swap, the D3 dismissal, and the D6 reframing were
  all decided by measurement that contradicted the written plan. Re-measure before trusting a
  claim in `README_plan.md` — several of its original numbers did not reproduce.
- **Never `error()` inside `Threads.@threads`.** One bad 6×6 would kill an entire sweep. Record
  and fall back; report the failure rate.
- **No number enters a paper table without an interval.**
