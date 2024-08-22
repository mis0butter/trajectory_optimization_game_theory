using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

rng = MersenneTwister( 1 )


## ============================================ ##
# init params 

params, players, game = init_game(  ) 

# compute all possible Δv solutions 
players = players_states( params, game, players, rng ) 

# save player state in game 
push!( game.p1_state, players[1] ) 
push!( game.p2_state, players[2] )  

# next game steps  

game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 

# plotting stuff 

# k_replan step of the game 
k = 2 
fig = plot_p1_p2_traj( game, params, k )  


## ============================================ ##



function dist_norm( game, params ) 

    tt_replan = params.tt_replan 
    k_max     = game.k_replan[end] 

    p1_r_hist = [] 
    p2_r_hist = [] 
    for kk = 1 : k_max 
    
        p1_chosen = game.p1_state[ kk ].chosen 
        p2_chosen = game.p2_state[ kk ].chosen  
    
        p1_r = game.p1_state[ kk ].X[ p1_chosen ][ 1 : tt_replan , 1 : 3 ] 
        p2_r = game.p2_state[ kk ].X[ p2_chosen ][ 1 : tt_replan , 1 : 3 ] 
    
        push!( p1_r_hist, p1_r ) 
        push!( p2_r_hist, p2_r ) 
    end 
    
    p1_r_hist = mapreduce( permutedims, hcat, p1_r_hist )' 
    p2_r_hist = mapreduce( permutedims, hcat, p2_r_hist )' 
    
    p1_r_norm = [ norm( p1_r_hist[ii,:] ) for ii in 1 : size(p1_r_hist, 1) ] 
    p2_r_norm = [ norm( p2_r_hist[ii,:] ) for ii in 1 : size(p2_r_hist, 1) ] 
    
    return p1_r_norm, p2_r_norm 
end 
