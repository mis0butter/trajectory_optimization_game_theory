
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
end 

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
    k_replan            # replan index 
    rv_E                # current vector for evader 
    rv_P                # current vector for pursuer 
    t_ref_E 
    rv_ref_E            # reference vector for polygon from rv_ref_E 
    p1_state            # player 1 state 
    p2_state            # player 2 state 
    params              # game parameters 

end 

export game_struct 


