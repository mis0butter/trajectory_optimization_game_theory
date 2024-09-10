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

# get costs 
games_costs, games_costs_mean = stage_cost_games_fn( games_vec )


## ============================================ ## 
 
tt_hist, _, _, _ = p_rv_ref_hist( game, params ) 

    # figure 
    fig = Figure( size = (600, 600) ) 

    # axis 1 
    # title_string = string( "p1 strategy = ", params.strategy, ", p2 strategy = ", params.p2_strategy, "\n", "mean cumsum norm of U vectors \n", "p1 = ", sprintf_stats.p1_Unorm_mean_end, ", p2 = ", sprintf_stats.p2_Unorm_mean_end )  
    title_string = "test"

    ax1 = Axis( fig[1,1], xlabel = "time (s)", title = title_string ) 
    # p1_ax1 = lines!( ax1, tt[ 1 : end - 1 ], stats.p1_Unorm_sum_mean, color = :blue ) 
    # p2_ax1 = lines!( ax1, tt[ 1 : end - 1 ], stats.p2_Unorm_sum_mean, color = :red ) 
    # Legend( fig[1,2], [ p1_ax1, p2_ax1 ], ["p1", "p2"] ) 
    for ii in eachindex(games_vec)
        lines!( ax1, tt[1:end-1], games_costs[ii,:][:], color = :blue, alpha = 0.1 ) 
        # lines!( ax1, tt[1:end-1], stats.p2_Unorm_sum_all[ii,:][:], color = :red,  alpha = 0.1 ) 
    end 



