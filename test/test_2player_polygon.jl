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

N_games  = 10 
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


## ============================================ ##

ii   = 9 
game = games_vec[ ii ] 
fig  = plot_game_stats( game, params ) 


## ============================================ ##

fig = plot_MC_stats( games_vec, params ) 


## ============================================ ##
# save 
 
data0  = [] .* zeros(1,4) 
header = [ "tt_hist", "p1_rv_hist", "p2_rv_hist", "rv_ref_hist" ] 
data_frame = DataFrame( data0, header ) 

for ii in eachindex(games_vec) 

    game = games_vec[ ii ] 
    
    tt_hist, p1_rv_hist, p2_rv_hist, rv_ref_hist = p_rv_ref_hist( game, params ) 

    data  = [ tt_hist , p1_rv_hist , p2_rv_hist , rv_ref_hist ] 
    push!( data_frame, data )
    
end 

save_folder   = string( "test/results/", params.strategy, "/")  
filename      = string( "games_N=", length(games_vec), ".csv" ) 
full_filename = string( save_folder, filename ) 

CSV.write( full_filename, data_frame, header=header ) 


## ============================================ ##

# reading this file is so stupid 

# load CSV file as data frame 
df = CSV.read( full_filename, DataFrame ) 

kk = 1 

str_chop = chop( df.tt_hist[kk]; head=1, tail=1 ) 
str_vec  = split( str_chop, ";") 

