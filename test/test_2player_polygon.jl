using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

rng = MersenneTwister( 1 )


## ============================================ ##
# init params 

params, players, game = init_game( rng ) 

## ============================================ ##
# next game steps  

game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 
game = prop_game_step( game, params, rng ) 

# plotting stuff 

# k_replan step of the game 
k = 1 
fig = plot_p1_p2_traj( game, params, k )  


## ============================================ ##

using GLMakie 

r_norm = dist_norm( game, params ) 
p1_U_norm, p2_U_norm = U_norm( game, params ) 
p1_Unorm_sum = cumsum( p1_U_norm ) 
p2_Unorm_sum = cumsum( p2_U_norm ) 

fig = Figure( resolution = (600, 800) )

ax1 = Axis( fig[1,1], xlabel = "time", title = "norm of U vectors" ) 
p1_ax1 = lines!( ax1, 1 : length(p1_U_norm), p1_U_norm, color = :blue ) 
p2_ax2 = lines!( ax1, 1 : length(p2_U_norm), p2_U_norm, color = :red ) 
Legend( fig[1,2], [ p1_ax1, p2_ax2 ], ["p1", "p2"] ) 

ax2 = Axis( fig[2,1], xlabel = "time", title = "sum of norm of U vectors" ) 
lines!( ax2, 1 : length(p1_Unorm_sum), p1_Unorm_sum, color = :blue ) 
lines!( ax2, 1 : length(p2_Unorm_sum), p2_Unorm_sum, color = :red ) 

ax3 = Axis( fig[3,1], xlabel = "time", title = "distance" ) 
lines!( ax3, 1 : length(r_norm), r_norm, color = :green )  

fig 

