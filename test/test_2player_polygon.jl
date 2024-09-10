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


    # figure 
    fig = Figure( size = (600, 600) ) 

    i_fig = 1 

    function plot_games_costs( fig, i_fig, games_vec ) 

        # get costs 
        games_costs, games_costs_mean, games_costs_std = stage_cost_games_fn( games_vec )

        # mean and mean-std strings 
        string_mean = @sprintf "%.3g" mean(games_costs_mean) 
        string_std_mean  = @sprintf "%.3g" mean(games_costs_std) 

        # get time vector 
        tt_hist, _, _, _ = p_rv_ref_hist( game ) 

        # mean +/- std 
        y_upper = games_costs_mean .+ games_costs_std 
        y_lower = games_costs_mean .- games_costs_std 

        # axis title 
        title_string = "stage costs (game value)"
        title_string = string( "mean cost = ", string_mean, ", mean std = ", string_std_mean )  

        # create axis 
        ax = Axis( fig[i_fig,1], xlabel = "time (s)", title = title_string )

        # plot mean 
        lines!( ax, tt_hist[1:end-1], games_costs_mean, color = :green )   

        # plot mean +/- std ribbons 
        fill_between!(ax, tt_hist[1:end-1], y_lower, y_upper, color = :green, alpha = 0.25 ) 

        # plot individual games 
        for ii in eachindex(games_vec)
            lines!( ax, tt_hist[1:end-1], games_costs[ii,:][:], color = :green, alpha = 0.1 ) 
        end 

        return fig 
    end 

fig 


    # # axis 2 
    # title_string = string( "mean player distance = ", sprintf_stats.dist_norm_mean_mean, ", mean std = ", @sprintf "%.3g" mean(stats.dist_norm_std) )
    # ax2 = Axis( fig[2,1], xlabel = "time (s)", title = title_string )  
    # y_upper = stats.dist_norm_mean .+ stats.dist_norm_std 
    # y_lower = stats.dist_norm_mean .- stats.dist_norm_std 
    # fill_between!(ax2, tt, y_lower, y_upper, color = :green, alpha = 0.25 )
    # lines!( ax2, tt, stats.dist_norm_mean, color = :green ) 
    # for ii in eachindex(games_vec)
    #     lines!( ax2, tt, stats.dist_rnorm_all[ii,:][:], color = :green, alpha = 0.1 ) 
    # end 


