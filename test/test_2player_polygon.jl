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

function U_norm( game ) 

    k_max = game.k_replan[end] 

    p1_U_hist = [] 
    p2_U_hist = [] 
    for kk = 1 : k_max 
    
        p1_chosen = game.p1_state[kk].chosen 
        p2_chosen = game.p2_state[kk].chosen  
    
        p1_U = game.p1_state[kk].U[p1_chosen] 
        p2_U = game.p2_state[kk].U[p2_chosen] 
    
        push!( p1_U_hist, p1_U ) 
        push!( p2_U_hist, p2_U ) 
    end 
    
    p1_U_hist = mapreduce( permutedims, hcat, p1_U_hist )' 
    p2_U_hist = mapreduce( permutedims, hcat, p2_U_hist )' 
    
    p1_U_norm = [ norm( p1_U_hist[ii,:] ) for ii in 1 : size(p1_U_hist, 1) ] 
    p2_U_norm = [ norm( p2_U_hist[ii,:] ) for ii in 1 : size(p2_U_hist, 1) ] 
    
    return p1_U_norm, p2_U_norm 
end 

