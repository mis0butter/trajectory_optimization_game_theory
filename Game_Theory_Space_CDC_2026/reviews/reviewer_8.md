Reviewer 8 of CDC 2026 submission 2010

# Comments to the author

The paper studies how learning an opponent’s behavior can
improve performance in an adversarial orbital setting. The
topic is interesting, and the simulation results are
useful. However, I found some parts of the mathematical
formulation unclear, especially the game-theoretic setup.
There are also several notation issues, some of which I
list below.

- The payoff matrix appears to be defined from the evader’s
perspective. If so, the evader and pursuer should solve
different optimization problems. The authors should write
the LPs for both players, or otherwise clearly explain how
each mixed strategy is computed using Eq(10).
- Eq(7) seems to compare the full state x_i(t_L) in R^6
with the position vertex r_vj, while they do not have the
same dimension.
- The indexing of the target vertices is confusing. In the
definition of A_jl, v_j seems to be the evader’s target
vertex and v_l the pursuer’s target vertex, but Eq(7)
uses r_vj as the terminal target for a generic agent i.
This indexing should be revised or explained more clearly.
- In Eq(14) and Eq(15), the same symbol is used for beliefs
over different opponent action spaces. It would be clearer
to use different symbols for these belief vectors.
- The term “dominant strategy” should be reconsidered. In
Table I, FP-greedy is not dominant for the evader, as
FP-mixed performs slightly better against an FP-greedy
pursuer. The authors should clarify this game-theoretic
term or use a different phrase.
- The game step notation (S) appears inconsistent across
the text, figures, and captions, creating confusion...
- Some references are incomplete or inconsistent in
formatting and should be revised.