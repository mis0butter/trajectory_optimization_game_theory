using trajectory_optimization_game_theory 

using LinearAlgebra: norm, dot 
using Statistics: mean, var, std 

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

rng = MersenneTwister( 1 )


## ============================================ ##
# init params 

params, players, game = init_game( rng ) ; 

## ============================================ ##
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

gameS = [] 
N_games = 10 
for jj = 1 : N_games 

    params, players, game = init_game( rng ) 

    N_replan = 10 
    for ii = 1 : N_replan - 1 
        print( "game: ", jj, " step: ", ii, "\n" ) 
        game = prop_game_step( game, params, rng ) 
    end 

    push!( gameS, game ) 

end 

## ============================================ ##

using GLMakie 

game = gameS[end-2] 

r_norm = dist_norm( game, params ) 
p1_U_norm, p2_U_norm = U_norm( game, params ) 
p1_Unorm_sum = cumsum( p1_U_norm ) 
p2_Unorm_sum = cumsum( p2_U_norm ) 

fig = Figure( resolution = (600, 800) )

ax1 = Axis( fig[1,1], xlabel = "time", title = "norm of U vectors" ) 
p1_ax1 = lines!( ax1, 1 : length(p1_U_norm), p1_U_norm, color = :blue ) 
p2_ax1 = lines!( ax1, 1 : length(p2_U_norm), p2_U_norm, color = :red ) 
Legend( fig[1,2], [ p1_ax1, p2_ax1 ], ["p1", "p2"] ) 

ax2 = Axis( fig[2,1], xlabel = "time", title = "cumsum of norm of U vectors" ) 
lines!( ax2, 1 : length(p1_Unorm_sum), p1_Unorm_sum, color = :blue ) 
lines!( ax2, 1 : length(p2_Unorm_sum), p2_Unorm_sum, color = :red ) 

ax3 = Axis( fig[3,1], xlabel = "time", title = "distance" ) 
lines!( ax3, 1 : length(r_norm), r_norm, color = :green )  

fig 

## ============================================ ##

dist_rnorm_all   = [] 
p1_Unorm_sum_all = [] 
p2_Unorm_sum_all = [] 
for ii in eachindex(gameS)

    game = gameS[ii] 

    # compute norm of U and distance vectors 
    dist_rnorm = dist_norm( game, params ) 
    p1_U_norm, p2_U_norm = U_norm( game, params ) 
    p1_Unorm_sum = cumsum( p1_U_norm ) 
    p2_Unorm_sum = cumsum( p2_U_norm ) 

    push!( dist_rnorm_all, dist_rnorm ) 
    push!( p1_Unorm_sum_all, p1_Unorm_sum ) 
    push!( p2_Unorm_sum_all, p2_Unorm_sum ) 

end 

dist_rnorm_all   = vv2m( dist_rnorm_all ) 
p1_Unorm_sum_all = vv2m( p1_Unorm_sum_all ) 
p2_Unorm_sum_all = vv2m( p2_Unorm_sum_all ) 

dist_norm_mean    = mean( dist_rnorm_all, dims = 1 )[:] 
p1_Unorm_sum_mean = mean( p1_Unorm_sum_all, dims = 1 )[:] 
p2_Unorm_sum_mean = mean( p2_Unorm_sum_all, dims = 1 )[:] 

T = params.tof * params.k_tt_replan 

# plot 
fig = Figure( resolution = (600, 600) ) 

title_string = string("games = ", length(gameS), "\n", "mean cumsum norm of U vectors")
ax1 = Axis( fig[1,1], xlabel = "time", title = title_string ) 
tt = ( 0 : length(p1_Unorm_sum_mean)-1 ) * T / length(p1_Unorm_sum_mean)  
p1_ax1 = lines!( ax1, tt, p1_Unorm_sum_mean, color = :blue ) 
p2_ax1 = lines!( ax1, tt, p2_Unorm_sum_mean, color = :red ) 
Legend( fig[1,2], [ p1_ax1, p2_ax2 ], ["p1", "p2"] ) 
for ii in eachindex(gameS)
    lines!( ax1, tt, p1_Unorm_sum_all[ii,:][:], color = :blue, alpha = 0.1 ) 
    lines!( ax1, tt, p2_Unorm_sum_all[ii,:][:], color = :red, alpha = 0.1 ) 
end 

ax2 = Axis( fig[2,1], xlabel = "time", title = "mean distance" ) 

tt = ( 0 : length(dist_norm_mean)-1 ) * T / length(dist_norm_mean)  
lines!( ax2, tt, dist_norm_mean, color = :green ) 
for ii in eachindex(gameS)
    lines!( ax2, tt, dist_rnorm_all[ii,:][:], color = :green, alpha = 0.1 ) 
end 

fig 
