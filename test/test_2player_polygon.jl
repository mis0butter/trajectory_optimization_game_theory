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
games_vec = run_MC_games(   rng, N_games, N_replan, "random",    "pure"    ) 
games_vec = run_MC_games(   rng, N_games, N_replan, "random",    "mixed"   ) 
games_vec = run_MC_games(   rng, N_games, N_replan, "mixed",     "pure"    ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "pure" ) 
# games_vec = run_MC_games( rng, N_games, N_replan, "random" ) 

fig = plot_MC_stats( games_vec ) 


## ============================================ ##
# load and plot games_vec 

games_vec = load_games_vec( 100, "pure", "pure" ) 
fig = plot_MC_stats( games_vec ) 


## ============================================ ## 

# get costs 
games_costs, games_costs_mean, games_costs_std = stage_cost_games_fn( games_vec )

tt_hist, _, _, _ = p_rv_ref_hist( game, params ) 

    # figure 
    fig = Figure( size = (600, 600) ) 

    # axis 1 
    title_string = "stage costs (game value)"

    ax1 = Axis( fig[1,1], xlabel = "time (s)", title = title_string ) 
    for ii in eachindex(games_vec)
        lines!( ax1, tt_hist[1:end-1], games_costs[ii,:][:], color = :blue, alpha = 0.1 ) 
    end 

fig 


