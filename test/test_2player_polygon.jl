using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean, var, std 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

rng = MersenneTwister( 1 )

using Infiltrator 


## ============================================ ##
# init params 

params, players, game = init_game( rng ) ; 

# ----------------------- #
# next game steps  

N_replan = 9  
for ii = 1 : N_replan - 1 
    print( "step: ", ii, "\n" ) 
    game = prop_game_step( game, params, rng ) 
end 

# ----------------------- #
# plotting stuff 

# k = 1 
fig = plot_p1_p2_traj( game, params, N_replan )  


## ============================================ ##
## ============================================ ##
# run multiple games 

N_games  = 10
N_replan = 10 

gameS = [] 
for jj = 1 : N_games 

    params, players, game = init_game( rng ) 

    for ii = 1 : N_replan - 1 
        print( "game: ", jj, " step: ", ii, "\n" ) 
        game = prop_game_step( game, params, rng ) 
    end 

    push!( gameS, game ) 

end 

## ============================================ ##

ii = 9 
fig = plot_game_stats( gameS, ii, params ) 

## ============================================ ##

fig = plot_MC_stats( gameS, params ) 

## ============================================ ##



