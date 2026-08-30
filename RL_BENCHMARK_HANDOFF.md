# Optimal-Control Benchmark Spec: N-Segment Impulsive Pursuit
### Handoff for RL implementation in a separate repo

Source of truth: `test/opt_Nseg_2player.jl` in the Julia repo
`trajectory_optimization_game_theory` (author: Junette Hsin). Everything below is
transcribed from that script and the functions it calls, so an RL agent can
reproduce the *exact same problem* and be compared against the optimal-control
solution head-to-head.

**Units everywhere: km, km/s, seconds, radians.**

---

## 1. Scenario in one paragraph

A **pursuer** and an **evader** start on two different LEO orbits, ~14,072 km
apart. The pursuer has a fixed time of flight `tof = 2000 s` to reach the
evader's position. The pursuer's control is discretized into **N = 30 impulsive
Δv burns**, one at the start of each of 30 equal-length coast segments
(`tof_N = tof/N = 66.667 s` each). In this particular script the **evader is
uncontrolled** — it coasts ballistically, so its terminal state is known in
advance and the problem reduces to a fixed-time, fixed-endpoint, impulsive
trajectory optimization. (The two-sided game version lives in
`src/game_theory/matrix_game_solver.jl`, which calls the same optimizer inside a
receding-horizon matrix game.)

---

## 2. Physical constants

| Symbol | Value | Meaning |
|---|---|---|
| `mu` | `398600.4415` | Earth gravitational parameter, km³/s² |
| `r`  | `6378.0` | Earth radius, km (used only to build the initial orbits) |
| `tof` | `2000` | total time of flight, s |
| `N` | `30` | number of segments / impulses |
| `tof_N` | `tof/N = 66.6667` | coast duration per segment, s |
| `Δv_max` | `2.0` | max magnitude of a *single* segment's Δv, km/s |

Dynamics are **two-body Keplerian only** — no J2, no drag, no third body.

---

## 3. Initial conditions

Classical orbital elements `[a, e, i, Ω, ω, ν]` (a in km, angles in radians):

```
# pursuer
kep0_P = [ 6778.0,  0.1, -20.0°,  10.0°,  20.0°,  30.0° ]
# evader
kep0_E = [ 6828.0,  0.2,  10.6°,  40.0°,   0.0°, 180.0° ]
```

Converted with `kep2cart` (perifocal → inertial via `R3(-Ω) * R1(-i) * R3(-ω)`,
i.e. the standard 3-1-3 PQW→IJK rotation), these are the exact Cartesian states
your RL environment should start from — **hard-code these to avoid any
convention mismatch**:

```
rv_0_P = [ 3137.24728672698,   5067.1066348626255, -1617.9745804804681,
             -7.000588510157188,  4.183857629142513,  -1.942121463958331 ]

rv_0_E = [-6276.661749139659,  -5266.744558727627,      0.0,
              3.9415673564967255, -4.697377057557039, -1.1475707834864168 ]
```

Derived quantities (sanity checks):

| Quantity | Pursuer | Evader |
|---|---|---|
| `‖r₀‖` (km) | 6175.41 | 8193.60 |
| `‖v₀‖` (km/s) | 8.3836 | 6.2384 |
| orbital period (s) | 5553.46 | 5615.02 |

Initial relative distance: **14,072.24 km**. `tof = 2000 s` ≈ **0.36 orbits** —
this is an aggressive intercept, which is why the optimal Δv comes out
physically huge (see §7). Keep that in mind before "fixing" anything.

### Target state

```julia
t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0)   # numerical 2-body ODE, saveat = 1 s
rv_f = rv_E[end, :]                                  # evader state at t = tof
```

`propagate_2Body` integrates `eom_2Body!` with DifferentialEquations.jl over
`tspan = (0, 2000)` saving every 1 s. Resulting target (also safe to hard-code):

```
rv_f = [ 5147.151264674194, -3251.8067208970774, -1085.3569414513438,
            3.237607251966227,   7.7156244491080175,  0.716655988454322 ]
```

Note: the *evader* is propagated with a **numerical ODE integrator**, while the
*pursuer* inside the optimizer is propagated with an **analytic Kepler
propagator** (§4). Small inconsistency in the original code; both are pure
two-body so the difference is integrator error only, but replicate it if you
want bit-comparable results.

---

## 4. Pursuer dynamics inside the optimizer

`prop_kepler_tof_Nseg(rv_0, Δv_vec, N, tof_N, mu)` — for `i = 1 … N`:

1. **Apply impulse**: `rv⁺ = [r ; v + Δv_i]` (instantaneous, position unchanged).
2. **Coast**: analytic Kepler propagation for `tof_N` seconds
   (`prop_kepler_tof`: cart→elements, advance true anomaly by solving Kepler's
   equation for the elapsed time, elements→cart; handles `e < 1` and `e ≥ 1`).

So the trajectory is: burn, coast 66.67 s, burn, coast, … 30 times. The final
state after the 30th coast is what gets compared to `rv_f`.

Δv vectors are expressed in the **inertial (ECI) frame**, not LVLH/RTN.

---

## 5. Decision variable

```
x = flatten(Δv_vec)   where Δv_vec is [N, 3] = [30, 3]   →   x ∈ ℝ⁹⁰
```

Row `i` is the 3-component Δv applied at the start of segment `i`.
`reshape(x, N, 3)` / `reshape(Δv, N*3, 1)` are used interchangeably throughout.

---

## 6. The cost function (this is the important part)

The script actually runs **`min_Δv_dist`** (`src/Opt/min_fns.jl:40`). Its
objective is an **unconstrained penalty form** — miss distance is folded
directly into the cost rather than imposed as a constraint:

```
J(x) =  Σ_{i=1}^{N} ‖Δv_i‖²        # "sum_norm_Δv"  — despite the name, sums SQUARES
      + ‖ r_final(x) − r_target ‖   # "miss_distance_prop_kepler_Nseg" — position only
```

Two details that are easy to get wrong:

- **`sum_norm_Δv` sums squared magnitudes**, not magnitudes:
  `sum_norm_Δv(x, N) = Σ_i Σ_j Δv[i,j]²`. The docstring says "sum of norms"; the
  code (`src/Opt/obj_fns.jl:160`) does `sum(Δv_vec[i,:].^2)`.
- **Miss distance uses position only** — `norm(rv_prop[1:3] − rv_f[1:3])`.
  Terminal *velocity* is never matched. This is an intercept, not a rendezvous.

The two terms are **added with weight 1.0 each and have inconsistent units**
(km²/s² + km). Since miss distance starts in the thousands of km, the optimizer
is effectively miss-distance-dominated early and fuel-shaped only near the end.
If you want a cleaner RL reward, keep this exact form for the benchmark
comparison and add a re-weighted variant as a separate experiment.

### Constraint

```
h_i(x) = ‖Δv_i‖ − Δv_max ≤ 0,   i = 1…N,   Δv_max = 2.0 km/s
```

(`constrain_Δv`, `src/Opt/obj_fns.jl:178`.) Per-segment, not cumulative. In the
converged solution the largest single burn is 0.59 km/s, so **this constraint is
inactive at the optimum** — in RL you can implement it as a simple action clip.

### The two sibling objectives (not used by default, both commented out in the script)

| Function | Objective | Constraints |
|---|---|---|
| `min_Δv` | `Σ‖Δv_i‖²` only | **equality**: miss distance = 0, plus the Δv_max inequality |
| `min_Δv_dist` ← **used** | `Σ‖Δv_i‖² + miss` | Δv_max inequality only |
| `max_Δv_dist` | `−Σ‖Δv_i‖² − miss` (evader-style) | Δv_max inequality; marked "seems to be not working" in the source |

---

## 7. Solver and baseline results (what RL must beat)

**Initial guess** (`lambert_init_guess`, `src/Lambert/Lambert.jl:75`): solve a
single-impulse Lambert problem (Battin's method) from `rv_0` to `rv_f` over
`tof`, take `Δv_lambert = v_lambert − v_0`, then split it evenly:
`Δv_i = Δv_lambert / N` for all i.

**Optimizer**: augmented Lagrangian outer loop (`min_aug_L` →
`min_aug_L_ineq`, `src/Opt/aug_L.jl`):
- λ initialized to 0, penalty p initialized to 10 for each of the N constraints, γ = 2
- inner unconstrained minimization by `Optim.optimize(fn, ∇fn, x, NelderMead())`
  with a `ForwardDiff` gradient supplied
- outer convergence when `‖x_{k+1} − x_k‖ < 1e-6` **and** all constraints strictly satisfied
- λ_i ← max(λ_i + p_i·h_i, 0), and p_i ← 2·p_i whenever h_i > 0

**Measured baseline** (this exact configuration, run 2026-08-24):

| Metric | Lambert initial guess | Optimized (`min_Δv_dist`) |
|---|---|---|
| Σ‖Δv_i‖ (km/s) | 14.860 | **10.644** |
| Σ‖Δv_i‖² (km²/s²) — the actual cost term | 7.361 | **4.160** |
| Terminal position miss (km) | 11,253.60 | **0.0955** |
| Max single-segment ‖Δv_i‖ (km/s) | — | 0.591 |
| Wall-clock solve time | — | ~2.0 s |

**Target for the RL agent: miss distance ≲ 0.1 km with Σ‖Δv_i‖² ≲ 4.16.**
Report both numbers — a policy can trivially cut Δv by missing, or cut miss
distance by burning more.

---

## 8. Suggested RL formulation (mapping, not prescription)

- **Episode**: exactly 30 steps, one per segment. Fixed horizon, fixed step time.
- **Action**: `a ∈ ℝ³` = Δv in the ECI frame, clipped to `‖a‖ ≤ 2.0` km/s.
  Consider normalizing to `‖a‖ ≤ 0.6` in practice — the optimum never exceeds
  that — but report results under the true 2.0 limit.
- **Observation**: at minimum the pursuer state `rv_P ∈ ℝ⁶`, the evader state
  `rv_E ∈ ℝ⁶` (or the relative state `rv_E − rv_P`), and remaining segments
  `(N − k)/N`. Normalize positions by ~7000 km and velocities by ~8 km/s.
- **Transition**: apply impulse to velocity, then Kepler-coast 66.667 s. Any
  accurate two-body propagator works; a universal-variable or Kepler-equation
  propagator matches the reference exactly.
- **Reward (to match the OC cost)**: per-step `−‖Δv_k‖²`, plus a terminal
  `−‖r_N − r_target‖`. Summed over an episode this is exactly `−J(x)`. A shaped
  variant (per-step reward on decreasing range) will train far better, but
  **evaluate with the unshaped `J`** so the numbers are comparable.
- **Termination**: only at k = 30 (no early intercept condition in the OC
  problem). If you add an intercept radius, note that as a deviation.

---

## 9. Deviations to watch for (things that will silently break comparability)

1. **Terminal velocity is not matched** — do not add a velocity term to the reward.
2. **Miss distance is position-only and unweighted** relative to the Δv² term.
3. The cost uses **squared** Δv magnitudes, so it is not a fuel-optimal
   (L1/Σ‖Δv‖) problem. It is closer to an energy-optimal formulation.
4. `kep2cart` uses the `R3(-Ω)·R1(-i)·R3(-ω)` convention and the pursuer's
   inclination is **negative** (−20°). Verify against the hard-coded Cartesian
   vectors in §3 before trusting your own converter.
5. `tof` is given as the integer `2000`; `tof_N = 2000/30 = 66.666…` is not a
   round number.
6. The evader target comes from a **numerical** integrator while the pursuer uses
   an **analytic** Kepler propagator.

---

## 10. Relevant source files (Julia repo)

| Path | Contents |
|---|---|
| `test/opt_Nseg_2player.jl` | the driver script this spec describes |
| `src/Opt/min_fns.jl` | `min_Δv`, `min_Δv_dist`, `max_Δv_dist`, `min_optim`, `min_bfgs` |
| `src/Opt/obj_fns.jl` | `sum_norm_Δv`, `miss_distance_prop_kepler_Nseg`, `constrain_Δv` |
| `src/Opt/aug_L.jl` | augmented Lagrangian machinery (`min_aug_L` and variants) |
| `src/Dyn/kepler.jl` | `prop_kepler_tof`, `prop_kepler_tof_Nseg` |
| `src/Dyn/propagator.jl` | `propagate_2Body`, `eom_2Body!`, `kep2cart`, `cart2kep` |
| `src/Dyn/prop_delta_v.jl` | `apply_Δv` |
| `src/Lambert/Lambert.jl` | `lambertbattin`, `prop_lambert_soln`, `lambert_init_guess` |
| `src/game_theory/matrix_game_solver.jl` | receding-horizon two-player game built on `min_Δv_dist` |
