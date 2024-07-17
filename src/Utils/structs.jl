
mutable struct player_struct 
    X 
    U 
    t 
    cost 
    weights 
    chosen 
end 

export player 

## ============================================ ##

mutable struct test_struct 
    A 
    B 
    C 
end 

export test_struct 

## ============================================ ##

mutable struct game_struct 

    tt          # timetag of the game 
    p1_state 
    p2_state 

end 

export game_struct 


