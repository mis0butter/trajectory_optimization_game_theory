using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean, var, std 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

using CSV, DataFrames 
using Infiltrator 

rng = MersenneTwister( 1 ) 


## ============================================ ##
# run single game 

rng = MersenneTwister( 1 ) 

# run_game( rng, N_replan = 10, p1_strategy = "mixed", p2_strategy = "mixed" ) 

N_replan = 10  
p1_strategy  = "mixed" 
p2_strategy  = "mixed" 
game, params = run_game( rng, N_replan, p1_strategy, p2_strategy ) 

# ----------------------- # 
# plotting stuff 

# k = 1 
fig = plot_p1_p2_traj( game, params, N_replan )  
fig = plot_game_stats( game, params ) 


## ============================================ ## 
## ============================================ ## 
# test running multiple games 

rng = MersenneTwister( 1 ) 

N_games  = 100 
N_replan = 10 

# p1_strategy  = "pure" 
# p2_strategy  = "mixed" 

# games_vec = run_MC_games( rng, N_games, N_replan, p1_strategy, p2_strategy ) 
# games_vec = run_MC_games(   rng, N_games, N_replan, "random",    "pure"    ) 
# games_vec = run_MC_games(   rng, N_games, N_replan, "random",    "mixed"   ) 
# games_vec = run_MC_games(   rng, N_games, N_replan, "mixed",     "pure"    ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "pure" ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "random" ) 

games_vec = run_MC_games_parallel( rng, N_games, N_replan, p1_strategy, p2_strategy ) 

fig = plot_MC_stats( games_vec ) 


# ## ============================================ ##
# # load and plot games_vec 

# games_vec = load_games_vec( 100, "pure", "pure" ) 
fig = plot_MC_stats( games_vec ) 


# ## ============================================ ## 


# figure 
fig = Figure( size = (600, 600) ) 

x_fig = 1 ; y_fig = 1 ; 

# ----------------------- #

# get stats 
stats = MC_stats( games_vec ) 
sprintf_stats = print_MC_stats( games_vec ) 

# temp strings for title string  
temp1 = @sprintf "%.3g" mean(stats.p1_ref_norm_std) 
temp2 = @sprintf "%.3g" mean(stats.p2_ref_norm_std) 

# axis title string 
title_string = string( 
    "mean player distance from reference orbit: ", 
    "\n p1 mean = ",  sprintf_stats.p1_ref_norm_mean_mean, 
    ", std = ", temp1, 
    "\n p2 mean = ",  sprintf_stats.p2_ref_norm_mean_mean, 
    ", std = ", temp2  
) 

# create axis 
ax3 = Axis( fig[x_fig, y_fig], xlabel = "time", title = title_string ) 

# plot reference norm mean 
p1_ax3 = lines!( ax3, tt, stats.p1_ref_norm_mean, color = :blue ) 
p2_ax3 = lines!( ax3, tt, stats.p2_ref_norm_mean, color = :red ) 

# plot individual games 
for ii in eachindex(games_vec)
    lines!( ax3, tt, stats.p1_ref_norm_all[ii,:][:], color = :blue, alpha = 0.1 ) 
    lines!( ax3, tt, stats.p2_ref_norm_all[ii,:][:], color = :red, alpha = 0.1 ) 
end 

# ----------------------- #

fig 



