# Comments to author (Associate Editor)

The paper studies a spacecraft pursuit–evasion problem using mixed trajectory strategies, fictitious play, and MPC-generated candidate trajectories within a non-cooperative orbital game framework. The reviewers acknowledged that 

- the application setting is interesting and that
- the simulation studies provide some useful practical insight into pursuit–evasion behavior.

The idea of combining trajectory generation with game-theoretic mixed strategies was viewed as potentially promising, and the simulation section was considered the strongest aspect of the paper.

However, the reviewers raised substantial concerns regarding the technical correctness, clarity, and theoretical grounding of the formulation. In particular, there appear to be

- inconsistencies in the payoff matrix definition and the associated linear program used to compute the mixed Nash strategy, potentially affecting the validity of the reported results.

Several key elements of the methodology are insufficiently defined or justified, including

- the interpretation and use of the Nash probability distribution,
- the meaning of the resulting mixed strategy when the MPC tracking is imperfect,
- and the rationale behind the greedy and fictitious-play-based sampling rules.

In addition, the paper lacks 

- theoretical analysis connecting the proposed discrete approximation to the original continuous control problem, and
- the claims regarding computational tractability are not adequately supported by solver details, runtime characterization, or guarantees.

While the overall direction is interesting, the current version requires substantial clarification and correction before the results can be considered reliable.