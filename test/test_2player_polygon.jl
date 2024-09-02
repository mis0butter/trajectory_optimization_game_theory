using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean, var, std 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

using CSV, DataFrames 
using Infiltrator 


## ============================================ ##
# init params 

rng = MersenneTwister( 1 )
params, players, game = init_game( rng, "mixed" ) ; 

# ----------------------- #
# next game steps  

N_replan = 10   
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

N_games  = 100 
N_replan = 10 

games_vec = [] 
for jj = 1 : N_games 

    params, players, game = init_game( rng, "random" ) 

    for ii = 1 : N_replan - 1 
        print( "game: ", jj, " step: ", ii, "\n" ) 
        game = prop_game_step( game, params, rng ) 
    end 

    push!( games_vec, game ) 

end 

fig = plot_MC_stats( games_vec, params ) 
print_MC_stats( games_vec, params ) 

save_games_vec_output( games_vec, params )  


## ============================================ ##
## ============================================ ##
# run multiple games 

N_games  = 100 
N_replan = 10 

games_vec = [] 
for jj = 1 : N_games 

    params, players, game = init_game( rng, "pure" ) 

    for ii = 1 : N_replan - 1 
        print( "game: ", jj, " step: ", ii, "\n" ) 
        game = prop_game_step( game, params, rng ) 
    end 

    push!( games_vec, game ) 

end 

fig = plot_MC_stats( games_vec, params ) 
print_MC_stats( games_vec, params ) 

save_games_vec_output( games_vec, params )  


## ============================================ ##
## ============================================ ##
# run multiple games 

N_games  = 100 
N_replan = 10 

games_vec = [] 
for jj = 1 : N_games 

    params, players, game = init_game( rng, "mixed" ) 

    for ii = 1 : N_replan - 1 
        print( "game: ", jj, " step: ", ii, "\n" ) 
        game = prop_game_step( game, params, rng ) 
    end 

    push!( games_vec, game ) 

end 

fig = plot_MC_stats( games_vec, params ) 
print_MC_stats( games_vec, params ) 

save_games_vec_output( games_vec, params )  







## ============================================ ##
## ============================================ ##
## ============================================ ##


ii   = 9 
game = games_vec[ ii ] 
fig  = plot_game_stats( game, params ) 


## ============================================ ##

fig = plot_MC_stats( games_vec, params ) 


## ============================================ ##
# save as CSV 

function save_games_vec_output( games_vec, params ) 

    # Define the folder and filename
    save_folder   = string("test/results/", params.strategy, "/")
    
    data0  = [] .* zeros(1,4) 
    header = [ "tt_hist", "p1_rv_hist", "p2_rv_hist", "rv_ref_hist" ] 
    data_frame = DataFrame( data0, header ) 

    for ii in eachindex(games_vec) 

        game_folder = string(save_folder, "game_", ii,"/") 

        if !isdir(game_folder)
            mkdir(game_folder)
        end 

        game = games_vec[ ii ] 
        
        tt_hist, p1_rv_hist, p2_rv_hist, rv_ref_hist = p_rv_ref_hist( game, params ) 

    
        CSV.write( string(game_folder, "tt_hist.csv"), DataFrame(tt_hist, :auto) ) 
        CSV.write( string(game_folder, "p1_rv_hist.csv"), DataFrame(p1_rv_hist, :auto) )  
        CSV.write( string(game_folder, "p2_rv_hist.csv"), DataFrame(p2_rv_hist, :auto) ) 
        CSV.write( string(game_folder, "rv_ref_hist.csv"), DataFrame(rv_ref_hist, :auto) ) 
        
    end 

end 


## ============================================ ##
# read CSV 

ii = 1 
game_folder = string(save_folder, "game_", ii,"/") 
filename    = string(game_folder, "tt_hist.csv")

tt_hist     = CSV.read( filename, DataFrame ) 




