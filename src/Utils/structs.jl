
mutable struct player_struct
    X           # contains candidate trajectories to each vertex
    U           # contains candidate control inputs to each vertex
    t           # contains candidate times to each vertex
    cost        # cost matrix
    weights     # weights for each vertex
    chosen      # chosen vertex
    rv_0_hist   # initial trajectory
    t_chosen    # chosen time
    rv_chosen   # chosen trajectory
    U_chosen    # chosen control input
    fp_belief   # fictitious play belief about opponent's vertex choices
    strategy_belief    # Bayesian belief about opponent's meta strategy
    tracked_strategies # List of tracked opponent meta strategies
    solve_info  # per-step solver diagnostics: timings, status, iterations
    learner_state      # persistent state for no-regret learners (regret/Hedge)
end

"""
    player_struct(; kwargs...)

Keyword constructor. Prefer this over the positional form.

Every field defaults to `[]` (or `nothing` for `learner_state`), so adding a
field here no longer requires editing every construction site in lockstep — which
is what previously made `IC.jl` and `play_games.jl` silently coupled.

`solve_info` is per-step and is NOT carried across game steps; it has the shape
`(traj = [...], lp = ...)`, where `traj` holds one record per candidate
trajectory optimization (12 per game step) and `lp` holds the single matrix-game
LP record. `learner_state` IS carried forward, like the belief fields.
"""
player_struct(;
    X                  = [],
    U                  = [],
    t                  = [],
    cost               = [],
    weights            = [],
    chosen             = [],
    rv_0_hist          = [],
    t_chosen           = [],
    rv_chosen          = [],
    U_chosen           = [],
    fp_belief          = [],
    strategy_belief    = [],
    tracked_strategies = [],
    solve_info         = (traj = [], lp = nothing),
    learner_state      = nothing,
) = player_struct(X, U, t, cost, weights, chosen, rv_0_hist, t_chosen,
                  rv_chosen, U_chosen, fp_belief, strategy_belief,
                  tracked_strategies, solve_info, learner_state)

export player_struct

## ====================================================================

mutable struct test_struct
    A
    B
    C
end

export test_struct

## ====================================================================

mutable struct game_struct

    tt                  # timetag of the game 
    i_game_step         # current game step index 
    rv_E                # current vector for evader 
    rv_P                # current vector for pursuer 
    t_ref_E
    rv_ref_E            # reference vector for polygon from rv_ref_E 
    p1_state            # player 1 state 
    p2_state            # player 2 state 
    params              # game parameters 

end

export game_struct


