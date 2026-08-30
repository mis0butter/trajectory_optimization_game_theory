Reviewer 7 of CDC 2026 submission 2010

# Comments to the author

The paper proposes an interesting framework for spacecraft
pursuit-evasion using candidate trajectory generation,
mixed strategies, and fictitious play. However, I have
several major concerns about the mathematical formulation
and the interpretation of the results.

Major comments:

1. Computational tractability is not sufficiently
justified. The paper claims that the proposed framework is
computationally tractable, but the main computational
burden is not the final 6 by 6 matrix game. The expensive
part is the generation of candidate trajectories through
the MPC problem. The paper does not specify enough
structure for the dynamics F, the feasible set, or the
optimization problem to determine whether this MPC problem
is convex, reliably solvable, or suitable for real-time
use. Solver details, runtime statistics, and assumptions on
F should be provided.
2. The cost/payoff convention is ambiguous and appears
inconsistent with the mixed-strategy LP. The continuous
game is written as a min-max problem, so the objective is
naturally a cost for the pursuer and a payoff for the
evader. The matrix A is then defined with evader strategies
as rows and pursuer strategies as columns. Therefore,
larger values of A should favor the evader and hurt the
pursuer. Under this convention, the evader is the row
maximizer and the pursuer is the column minimizer. However,
the LP in Eq. (10), with max 1^T z subject to A^T z <= 1,
appears to correspond to a row-player minimization
formulation, not the evader’s maximizing strategy. The
paper should clearly define whether A is an evader payoff
or pursuer cost, and should provide consistent primal/dual
LPs for both players.
3. The computation of mixed Nash strategies is
underexplained. The paper presents one LP and then states
that each agent samples from its respective weights.
However, in a zero-sum matrix game, the row and column
players generally have different mixed strategies. The
paper should clearly explain how the evader’s and pursuer’s
weights are separately computed, and why the presented LP
gives the desired strategy for each player.
4. The strategy comparison requires stronger justification.
The paper compares Greedy, Mixed, Random, FP-greedy, and
FP-mixed strategies, but some of these choices are not
sufficiently motivated. In particular, the Greedy strategy
selects the pure action with the largest Nash equilibrium
weight. This is not generally a rational use of a mixed
Nash equilibrium, because the equilibrium property depends
on randomizing according to the full distribution. If this
is intended as a heuristic baseline, the paper should state
this explicitly.
5. The FP-mixed strategy is underdefined. It appears to be
a heuristic stochastic version of fictitious play, but the
paper does not justify or specify it mathematically. The
paper says FP-mixed samples from c_exp, but c_exp is an
expected payoff/cost vector, not a probability
distribution. Its entries may be unnormalized or even
negative, especially for the pursuer where c_exp = -A^T
q_hat. The authors should explain how c_exp is converted
into valid sampling probabilities.
6. The term “game value” is used ambiguously. In zero-sum
game theory, game value usually refers to the minimax
equilibrium value. In the results section, however, game
value is described as the mean stage cost and is sometimes
averaged across opponents. This is an empirical performance
metric, not necessarily the game-theoretic value. The paper
should either reserve “game value” for the minimax value or
rename the reported metric as realized mean stage cost or
average realized payoff.

Minor comments:

1. The stage cost is defined differently in Eq. (2) and Eq.
(9). Eq. (2) uses only inter-player distance, while Eq. (9)
adds smoothing and relative control effort. The paper
should use distinct notation for physical distance and
augmented stage cost.
2. The empirical belief update should be clearly presented
as part of the fictitious-play strategies only. As written,
it may be read as applying to all strategies. The notation
in this subsection is also confusing, since the text
introduces n but updates f, and the indexing in the
empirical frequency definition should be clarified.
3. The reported metrics are partly overlapping. The mean
game value is based on the stage cost, which already
depends on distance and control effort. Reporting mean
distance and cumulative delta-V is useful, but the paper
should explain that these metrics partly decompose the same
underlying objective rather than treating them as fully
independent evidence.

###