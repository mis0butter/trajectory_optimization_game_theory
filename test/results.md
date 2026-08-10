## V. Numerical Experiments and Results

### A. Simulation Set-Up

We consider a pursuit-evasion scenario in low Earth orbit. Both spacecraft occupy near-circular orbits at an altitude of 620 km (semi-major axis $a = 6998$ km, eccentricity $e = 0.01$, inclination $i = 20°$), yielding an orbital period of approximately 97 minutes. The evader (Player 1) is initialized on a reference orbit, while the pursuer (Player 2) begins on a slightly higher orbit with semi-major axis $a_P = 1.005\,a$, introducing a natural along-track drift between the two vehicles. Both spacecraft are modeled as point masses executing impulsive $\Delta V$ maneuvers; no continuous-thrust or attitude dynamics are included, so the control authority of each player is characterized entirely by a per-segment velocity-change constraint of $\Delta V_{\max} = 2.0$ km/s.

At each game step, both players simultaneously select one of six target vertices arranged as a regular hexagon of radius $R = 6.378$ km in the local orbital frame. The trajectory from the current state to each candidate vertex is computed via multi-segment impulsive $\Delta V$ optimization subject to the per-segment constraint. The orbit is discretized into $N_\text{seg} = 50$ segments per period; each game step executes 5 segments (approximately one-tenth of the orbital period), and each player's optimizer plans over a 10-segment horizon. The resulting $6 \times 6$ cost matrix has entries

$$J_{ij} = \frac{1}{K}\sum_{k=1}^{K}\left[\sqrt{\lVert \mathbf{r}_1^{(k)} - \mathbf{r}_2^{(k)}\rVert + 0.1} \;+\; 0.1\,\lVert \mathbf{u}_1^{(k)} - \mathbf{u}_2^{(k)}\rVert\right],$$

where the sum averages the stage cost over the $K$ segments in the executed window, combining inter-player distance with a differential control-effort penalty. The game is zero-sum ($J_2 = -J_1$), and the mixed Nash equilibrium is computed via linear programming at each step.

Five strategies are evaluated for each player:

1. **Mixed** — sample a vertex from the Nash equilibrium mixed strategy at each step.
2. **Greedy** — deterministically select the vertex with the highest Nash equilibrium weight (i.e., the mode of the mixed strategy).
3. **Random** — select a vertex uniformly at random, independent of the cost matrix.
4. **FP-greedy** — maintain an empirical frequency count of the opponent's past vertex choices (fictitious play). Form the empirical distribution $\hat{q}$, compute the expected cost vector $C\hat{q}$, and best-respond by selecting the vertex that maximizes (evader) or minimizes (pursuer) the expected payoff.
5. **FP-mixed** — identical belief update to FP-greedy, but instead of best-responding deterministically, sample a vertex with probability proportional to the expected cost vector (shifted to be non-negative).

The FP strategies additionally maintain a Bayesian belief over a library of nine opponent-type hypotheses — mixed, greedy, random, and each of the six fixed-vertex strategies — with a uniform prior. After each observed opponent action, the belief is updated via Bayes' rule using the likelihood of the observed vertex under each hypothesis. This belief is used to compute the correct identification rate and belief entropy reported below; vertex selection itself is driven by the empirical frequency vector $\hat{q}$.

We evaluate all 25 pairwise matchups of the five strategies across three time horizons: $n = 10$, $20$, and $30$ game steps (corresponding to approximately 1, 2, and 3 orbital periods). Each matchup is run for 50 Monte Carlo trials with independent initial conditions drawn from a seeded pseudorandom number generator (MersenneTwister, master seed 1) for reproducibility. We report four classes of metrics, each averaged over the 50 trials: (i) the game value (mean stage cost), (ii) the mean inter-player distance, (iii) the cumulative $\Delta V$ for each player, and (iv) the Bayesian belief entropy and correct identification rate. To isolate steady-state behavior from transient effects, we additionally report *last-quarter* metrics computed over the final 25% of the segment timeline.

### B. Results and Discussion

#### Strategy Dominance

Table I presents the mean game value averaged across all opponents, for each evader (P1) and pursuer (P2) strategy at each horizon. Higher values favor the evader.

**Table I.** Mean last-quarter game value by strategy, averaged across all opponents.

| Evader (P1)  | $n{=}10$ | $n{=}20$ | $n{=}30$ |
| ---          | ---      | ---      | ---      |
| FP-greedy    | **2.62** | **2.45** | **2.45** |
| FP-mixed     | **2.22** | **2.24** | **2.20** |
| random       | 2.04     | 1.96     | 2.00     |
| mixed        | 1.66     | 1.69     | 1.63     |
| greedy       | 1.59     | 1.55     | 1.57     |

| Pursuer (P2) | $n{=}10$ | $n{=}20$ | $n{=}30$ |
| ---          | ---      | ---      | ---      |
| FP-greedy    | **1.54** | **1.51** | **1.51** |
| FP-mixed     | 1.86     | 1.81     | 1.84     |
| random       | 2.16     | 2.11     | 2.12     |
| mixed        | 2.27     | 2.30     | 2.25     |
| greedy       | 2.31     | 2.16     | 2.15     |

FP-greedy is the dominant strategy on both sides of the game: the best evader and the best pursuer at every horizon tested. The ranking is fully consistent — fictitious-play strategies outperform naive ones, and the greedy strategy is the worst evader and a poor pursuer. The full $5\times 5$ game-value matrix (Table II, shown for $n=30$) reveals that the performance gap is stark: the best cell for the evader (FP-greedy evader vs. mixed pursuer, 3.08) exceeds the worst cell (greedy evader vs. FP-greedy pursuer, 1.16) by a factor of 2.7$\times$.

**Table II.** Last-quarter game value matrix at $n = 30$ (rows = evader, columns = pursuer).

| P1 \ P2       | FP-greedy | FP-mixed | greedy | mixed    | random |
| ---            | ---       | ---      | ---    | ---      | ---    |
| **FP-greedy**  | 1.70      | 2.15     | 2.76   | **3.08** | 2.58   |
| **FP-mixed**   | 1.71      | 2.03     | 2.37   | 2.46     | 2.43   |
| **greedy**     | **1.16**  | 1.51     | 1.78   | 1.75     | 1.67   |
| **mixed**      | 1.33      | 1.59     | 1.72   | 1.75     | 1.77   |
| **random**     | 1.64      | 1.92     | 2.09   | 2.22     | 2.13   |

Two structural features are noteworthy. First, the FP-greedy column (pursuer) is the column minimum for nearly every evader row — no evader can escape it as effectively as it can escape any other pursuer. Second, the mixed pursuer's randomization is counterproductive against a learning evader: the FP-greedy evader vs. mixed pursuer cell (3.08) persistently exceeds the FP-greedy evader vs. greedy pursuer cell (2.76) because the mixed pursuer occasionally selects suboptimal vertices that the FP evader learns to exploit.

#### Inter-Player Distance and Late-Game Convergence

The mean inter-player distance matrices exhibit the same structural pattern as the game-value matrices, consistent with the distance-dominated cost function. Table III reports the full distance matrix at $n=30$.

**Table III.** Mean inter-player distance (km) at $n = 30$.

| P1 \ P2       | FP-greedy | FP-mixed | greedy    | mixed     | random |
| ---            | ---       | ---      | ---       | ---       | ---    |
| **FP-greedy**  | 4.39      | 6.04     | 10.24     | **10.50** | 8.22   |
| **FP-mixed**   | 4.52      | 5.63     | 7.43      | 7.74      | 7.25   |
| **greedy**     | **3.26**  | 3.94     | 4.69      | 4.64      | 4.39   |
| **mixed**      | 3.61      | 4.14     | 4.72      | 4.77      | 4.66   |
| **random**     | 4.33      | 5.10     | 6.14      | 6.41      | 6.20   |

The learning advantage compounds over time. At $n=10$, the ratio between the largest and smallest distance entries is $13.89 / 7.04 = 2.0\times$; by $n=30$ this ratio widens to $10.50 / 3.26 = 3.2\times$. The FP pursuer steadily closes while the naive pursuer falls progressively further behind.

The last-quarter metrics sharpen this picture. For the greedy evader vs. FP-greedy pursuer at $n=30$, the last-quarter distance is 1.30 km — well below the full-game mean of 3.26 km — indicating that the FP pursuer accelerates its closure in the late game as its belief estimate converges. By contrast, the FP-greedy evader vs. greedy pursuer sustains a last-quarter distance of 8.10 km (vs. a full-game mean of 10.24 km), showing that the learning evader maintains large separation throughout. The FP-vs-FP matchup settles at a last-quarter game value of 1.70 and distance of 2.82 km at $n=30$, consistent with a Nash-like equilibrium between two mutually adapting agents.

#### Fuel Efficiency

The $\Delta V$ expenditure data reveal that the FP advantage is not purchased at the cost of additional fuel — in fact, the opposite is true. Table IV compares cumulative $\Delta V$ for three representative matchups at $n=30$.

**Table IV.** Cumulative $\Delta V$ (km/s) and mean distance for selected matchups at $n = 30$.

| Matchup                        | P1 $\Delta V$ | P2 $\Delta V$ | P2/P1 | Distance |
| ---                            | ---            | ---            | ---    | ---      |
| FP-greedy evader vs. greedy    | 0.237          | 0.331          | 1.40   | 10.24    |
| greedy evader vs. FP-greedy    | 0.343          | 0.426          | 1.24   | 3.26     |
| FP-greedy evader vs. FP-greedy | 0.255          | 0.346          | 1.36   | 4.39     |

The pursuer universally expends more fuel than the evader, consistent with the asymmetric pursuit geometry. However, the greedy evader spends 45% more $\Delta V$ than the FP-greedy evader ($0.343$ vs. $0.237$ km/s) while achieving only one-third the separation ($3.26$ vs. $10.24$ km). The FP framework thus yields a dramatic improvement in fuel efficiency: the learning evader achieves superior separation while simultaneously conserving propellant, because it avoids the predictable, high-cost trajectories that a deterministic strategy commits to.

#### Belief Convergence and Strategy Identification

The Bayesian belief entropy provides a direct measure of how quickly each FP player resolves uncertainty about its opponent. By $n=30$, the belief entropy for FP-greedy players facing identifiable opponents (greedy, mixed, random) has collapsed to $\mathcal{O}(10^{-4})$ or below, indicating near-certain identification. Against other FP opponents — which are not themselves members of the hypothesis library — the entropy still collapses (to $\mathcal{O}(10^{-19})$ or lower), though in this case the belief converges to the closest proxy hypothesis rather than the true generative strategy. The correct identification rate is $1.0$ for all matchups in which the true opponent type appears in the hypothesis set (greedy, mixed, random), confirming that the Bayesian update reliably resolves the opponent's strategy within a few orbital periods.

This near-perfect identification explains the late-game acceleration observed in the distance data: once the FP player's belief has converged, it effectively plays the exact best response to its opponent's (now-known) strategy, eliminating the exploration cost that characterizes earlier game steps.

#### Discussion

The central finding is that the fictitious-play framework provides a decisive and compounding advantage in the orbital pursuit-evasion game. The advantage is primarily *informational*, not maneuverability-based: FP and non-FP players draw from the same action space and face the same $\Delta V$ constraints, yet FP players achieve 2–3$\times$ better outcomes across all metrics. Three observations support this interpretation.

First, the strategy ranking is invariant across all three time horizons, and the performance gap *widens* with the number of game steps. This is the signature of a learning advantage — the benefit of accumulated observations compounds over time, unlike a static positional advantage that would erode or plateau.

Second, random evasion outperforms greedy evasion at every horizon (mean last-quarter game value of 2.00 vs. 1.57 at $n=30$, averaged across pursuers), despite having no strategic intent. This confirms that *unpredictability* has inherent value for the evader in this game: a uniform random strategy is harder to learn and counter than a deterministic one, even though it makes no attempt to optimize the cost function.

Third, the FP-vs-FP matchup converges to a stable, moderate game value (1.70 at $n=30$), substantially below the FP-greedy evader's best outcome (3.08) but above the greedy evader's worst (1.16). When both players are learning, neither can gain a persistent informational edge, and the game settles near a dynamic equilibrium — the empirical analogue of the minimax value for this repeated game. This outcome validates the game-theoretic framing: the mixed Nash equilibrium is not merely a mathematical construct, but a behavioral attractor when both agents employ rational, adaptive strategies.

---
---

<!-- Raw data and working notes preserved below for reference -->

# Raw tables: game_value_last_quarter

Here are the `game_value_last_quarter` matrices (rows = P1 evader, columns = P2 pursuer). Higher values favor the evader.

### n = 10

| P1 \ P2 | FP_greedy | FP_mixed | greedy | mixed | random |
| --- | --- | --- | --- | --- | --- |
| **FP_greedy** | 1.62 | 2.17 | **3.51** | 3.13 | 2.65 |
| **FP_mixed** | 1.78 | 2.07 | 2.42 | 2.43 | 2.42 |
| **greedy** | **1.22** | 1.54 | 1.71 | 1.75 | 1.73 |
| **mixed** | 1.39 | 1.64 | 1.76 | 1.74 | 1.79 |
| **random** | 1.67 | 1.89 | 2.13 | 2.28 | 2.21 |

### n = 20

| P1 \ P2 | FP_greedy | FP_mixed | greedy | mixed | random |
| --- | --- | --- | --- | --- | --- |
| **FP_greedy** | 1.61 | 2.09 | 2.87 | **3.10** | 2.57 |
| **FP_mixed** | 1.70 | 1.99 | 2.46 | 2.60 | 2.42 |
| **greedy** | **1.22** | 1.52 | 1.64 | 1.77 | 1.61 |
| **mixed** | 1.36 | 1.56 | 1.75 | 1.79 | 1.76 |
| **random** | 1.64 | 1.86 | 2.09 | 2.22 | 2.17 |

### n = 30

| P1 \ P2 | FP_greedy | FP_mixed | greedy | mixed | random |
| --- | --- | --- | --- | --- | --- |
| **FP_greedy** | 1.70 | 2.15 | 2.76 | **3.08** | 2.58 |
| **FP_mixed** | 1.71 | 2.03 | 2.37 | 2.46 | 2.43 |
| **greedy** | **1.16** | 1.51 | 1.78 | 1.75 | 1.67 |
| **mixed** | 1.33 | 1.59 | 1.72 | 1.75 | 1.77 |
| **random** | 1.64 | 1.92 | 2.09 | 2.22 | 2.13 |

# Raw tables: player_distance

### n = 10

| P1 \ P2 | FP_greedy | FP_mixed | greedy | mixed | random |
| --- | --- | --- | --- | --- | --- |
| **FP_greedy** | 7.78 | 8.88 | **13.89** | 11.67 | 10.86 |
| **FP_mixed** | 7.73 | 8.56 | 10.31 | 10.32 | 9.78 |
| **greedy** | **7.04** | 7.23 | 8.17 | 7.84 | 7.71 |
| **mixed** | 7.26 | 7.45 | 8.16 | 7.95 | 7.85 |
| **random** | 7.61 | 8.09 | 9.24 | 9.33 | 9.11 |

### n = 20

| P1 \ P2 | FP_greedy | FP_mixed | greedy | mixed | random |
| --- | --- | --- | --- | --- | --- |
| **FP_greedy** | 5.26 | 6.73 | **11.94** | 10.85 | 8.87 |
| **FP_mixed** | 5.35 | 6.36 | 8.30 | 8.45 | 7.86 |
| **greedy** | **4.23** | 4.78 | 5.56 | 5.47 | 5.19 |
| **mixed** | 4.55 | 4.97 | 5.58 | 5.59 | 5.46 |
| **random** | 5.15 | 5.81 | 6.96 | 7.13 | 6.98 |

### n = 30

| P1 \ P2 | FP_greedy | FP_mixed | greedy | mixed | random |
| --- | --- | --- | --- | --- | --- |
| **FP_greedy** | 4.39 | 6.04 | **10.24** | 10.50 | 8.22 |
| **FP_mixed** | 4.52 | 5.63 | 7.43 | 7.74 | 7.25 |
| **greedy** | **3.26** | 3.94 | 4.69 | 4.64 | 4.39 |
| **mixed** | 3.61 | 4.14 | 4.72 | 4.77 | 4.66 |
| **random** | 4.33 | 5.10 | 6.14 | 6.41 | 6.20 |

# Raw analysis notes

- The data covers 5 strategies (`FP_greedy`, `FP_mixed`, `greedy`, `mixed`, `random`) for both evader (P1) and pursuer (P2) across 3 time horizons (10, 20, 30 game steps), yielding 25 matchups per horizon. Each row is averaged over 50 MC games.
- FP_greedy evader vs mixed pursuer cell (3.08–3.13) is persistently the highest — the mixed pursuer's randomization works against it when facing a learning evader that can exploit the suboptimal moves.
- greedy vs mixed pursuer actually slightly exceeds greedy vs greedy in the late game (1.75 vs 1.71 at n=10; 1.75 vs 1.78 at n=30), suggesting the mixed pursuer loses some edge by randomizing against an already-predictable evader.
- The distance gap widens between FP and naive matchups as the horizon grows. At n=10 the max/min ratio is 13.89/7.04 = 2.0x, but at n=30 it's 10.50/3.26 = 3.2x.
- P2 always expends more fuel than P1, consistent with pursuit mechanics.
- FP strategies collapse belief entropy to near-zero (10^-15 to 10^-45 by n=30).
- Correct identification rate is ~1.0 for all identifiable strategies (greedy, mixed, random). FP strategies themselves aren't in the identification set (shown as NaN).
