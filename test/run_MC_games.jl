using trajectory_optimization_game_theory

using LinearAlgebra: norm, dot
using Statistics: mean, var, std

using StatsBase: ProbabilityWeights, sample
using Random: MersenneTwister

using CSV, DataFrames
using Infiltrator


## ==================================================================== 
## run single game 
## ==================================================================== 

rng = MersenneTwister(1)

# run_game( rng, n_game_steps = 10, p1_strategy = "mixed", p2_strategy = "mixed" ) 

n_game_steps = 10
p1_strategy = 1
p2_strategy = "FP_mixed"
game, params = run_game(rng, n_game_steps, p1_strategy, p2_strategy)

# ---------------------------------- 
# plotting stuff 

# k = 1 
fig = plot_p1_p2_traj(game, params, n_game_steps)
fig_stats = plot_game_stats(game, params)


## ==================================================================== 
## test running multiple games 
## ==================================================================== 

rng = MersenneTwister(1)

N_games = 100
n_game_steps = 10

# greedy, mixed, random, FP_greedy, FP_mixed, Meta_greedy, Meta_mixed
p1_strategy = "greedy"
p2_strategy = "Meta_greedy"

# games_vec = run_MC_games( rng, N_games, n_game_steps, "random" ) 

games_vec = run_MC_games_parallel(rng, N_games, n_game_steps, p1_strategy, p2_strategy)

## ==================================================================== 
## run all cases 
## ==================================================================== 

rng = MersenneTwister(1)

N_games = 100
n_game_steps = 10

for p1_strategy in ["greedy", "mixed", "random", "FP_greedy", "FP_mixed", "Meta_greedy", "Meta_mixed"]
    for p2_strategy in ["greedy", "mixed", "random", "FP_greedy", "FP_mixed", "Meta_greedy", "Meta_mixed"]
        games_vec = run_MC_games_parallel(rng, N_games, n_game_steps, p1_strategy, p2_strategy)
    end
end



## ==================================================================== 
## plot stats 
## ==================================================================== 

fig = plot_MC_stats(games_vec)


## ==================================================================== 
## save as csvs 
## ==================================================================== 


"Extract rv histories from games_vec and save as CSV files"
function save_rv_histories_csv(games_vec, p1_strategy, p2_strategy, output_dir="test/game_rv_hists/")

    # Update output directory name based on strategies 
    output_dir = "$(output_dir)p1_$(p1_strategy)_p2_$(p2_strategy)"

    # Create output directory if it doesn't exist
    if !isdir(output_dir)
        mkdir(output_dir)
    end

    N_games = length(games_vec)

    # ---------------------------------- 
    # loop through each game and extract rv histories 
    # ---------------------------------- 

    for i_game = 1:N_games

        game = games_vec[i_game]

        p1_state = game.p1_state
        p2_state = game.p2_state

        # Initialize vectors for flattened data 
        t_vec = Float64[]
        p1_x = Float64[]
        p1_y = Float64[]
        p1_z = Float64[]
        p1_vx = Float64[]
        p1_vy = Float64[]
        p1_vz = Float64[]
        p1_ux = Float64[]
        p1_uy = Float64[]
        p1_uz = Float64[]
        p2_x = Float64[]
        p2_y = Float64[]
        p2_z = Float64[]
        p2_vx = Float64[]
        p2_vy = Float64[]
        p2_vz = Float64[]
        p2_ux = Float64[]
        p2_uy = Float64[]
        p2_uz = Float64[]

        # Track cumulative time offset across replan segments 
        t_offset = 0.0

        # ---------------------------------- 
        # loop through each replan segment 
        # ---------------------------------- 

        for i_replan = 1:length(p1_state)

            p1_t = p1_state[i_replan].t_chosen[1:end-1]
            p1_rv = p1_state[i_replan].rv_chosen[1:end-1, :]
            p1_u = p1_state[i_replan].U_chosen[1:end-1, :]

            p2_t = p2_state[i_replan].t_chosen[1:end-1]
            p2_rv = p2_state[i_replan].rv_chosen[1:end-1, :]
            p2_u = p2_state[i_replan].U_chosen[1:end-1, :]

            # ---------------------------------- 
            # Loop through each time point in this replan segment 
            # ---------------------------------- 
            # rv_chosen is N_pts × 6 matrix (rows = time pts, cols = x,y,z,vx,vy,vz)

            N_pts = length(p1_t)
            for i_pt = 1:N_pts
                push!(t_vec, p1_t[i_pt] + t_offset)
                push!(p1_x, p1_rv[i_pt, 1])
                push!(p1_y, p1_rv[i_pt, 2])
                push!(p1_z, p1_rv[i_pt, 3])
                push!(p1_vx, p1_rv[i_pt, 4])
                push!(p1_vy, p1_rv[i_pt, 5])
                push!(p1_vz, p1_rv[i_pt, 6])
                push!(p2_x, p2_rv[i_pt, 1])
                push!(p2_y, p2_rv[i_pt, 2])
                push!(p2_z, p2_rv[i_pt, 3])
                push!(p2_vx, p2_rv[i_pt, 4])
                push!(p2_vy, p2_rv[i_pt, 5])
                push!(p2_vz, p2_rv[i_pt, 6])
                push!(p1_ux, p1_u[i_pt, 1])
                push!(p1_uy, p1_u[i_pt, 2])
                push!(p1_uz, p1_u[i_pt, 3])
                push!(p2_ux, p2_u[i_pt, 1])
                push!(p2_uy, p2_u[i_pt, 2])
                push!(p2_uz, p2_u[i_pt, 3])
            end

            # Update offset for next segment (add the final time of this segment)
            t_offset += p1_state[i_replan].t_chosen[end]

        end

        # Create DataFrame with vector columns 
        df = DataFrame(
            t=t_vec,
            p1_x=p1_x, p1_y=p1_y, p1_z=p1_z,
            p1_vx=p1_vx, p1_vy=p1_vy, p1_vz=p1_vz,
            p1_ux=p1_ux, p1_uy=p1_uy, p1_uz=p1_uz,
            p2_x=p2_x, p2_y=p2_y, p2_z=p2_z,
            p2_vx=p2_vx, p2_vy=p2_vy, p2_vz=p2_vz,
            p2_ux=p2_ux, p2_uy=p2_uy, p2_uz=p2_uz,
        )
        CSV.write(joinpath(output_dir, "game_$(i_game).csv"), df)

        println("Saved game_", i_game, ".csv")

    end

    println("Saved ", N_games, " game files to ", output_dir)

end

export save_rv_histories_csv

# ---------------------------------- 

save_rv_histories_csv(games_vec, p1_strategy, p2_strategy)

## ==================================================================== 
## save reference orbit history as csv 
## ==================================================================== 

function save_ref_orbit_history_csv(game)

    t_ref_E_hist = game.t_ref_E
    rv_ref_E_hist = game.rv_ref_E

    t_vec = Float64[]
    ref_x = Float64[]
    ref_y = Float64[]
    ref_z = Float64[]

    for i_replan = 1:length(t_ref_E_hist)

        t_ref_E = t_ref_E_hist[i_replan]
        rv_ref_E = rv_ref_E_hist[i_replan]

        for i_pt = 1:5
            push!(t_vec, t_ref_E[i_pt])
            push!(ref_x, rv_ref_E[i_pt, 1])
            push!(ref_y, rv_ref_E[i_pt, 2])
            push!(ref_z, rv_ref_E[i_pt, 3])
        end
    end

    df = DataFrame(t=t_vec, ref_x=ref_x, ref_y=ref_y, ref_z=ref_z)

    CSV.write("test/game_rv_hists/ref_orbit_history.csv", df)
end

export save_ref_orbit_history_csv

# ---------------------------------- 

save_ref_orbit_history_csv(game)

## ==================================================================== 
## sanity check: plot trajectories from csv 
## ==================================================================== 

using Plots

"Plot p1 and p2 xyz trajectories from a game CSV file"
function plot_game_csv(csv_path)

    df = CSV.read(csv_path, DataFrame)

    # 3D trajectory plot 
    fig_3d = plot(df.p1_x, df.p1_y, df.p1_z,
        label="P1 (evader)", lw=2, color=:blue)
    plot!(fig_3d, df.p2_x, df.p2_y, df.p2_z,
        label="P2 (pursuer)", lw=2, color=:red)
    scatter!(fig_3d, [df.p1_x[1]], [df.p1_y[1]], [df.p1_z[1]],
        label="P1 start", ms=6, color=:blue)
    scatter!(fig_3d, [df.p2_x[1]], [df.p2_y[1]], [df.p2_z[1]],
        label="P2 start", ms=6, color=:red)
    scatter!(fig_3d, [df.p1_x[end]], [df.p1_y[end]], [df.p1_z[end]],
        label="P1 end", ms=6, marker=:star, color=:blue)
    scatter!(fig_3d, [df.p2_x[end]], [df.p2_y[end]], [df.p2_z[end]],
        label="P2 end", ms=6, marker=:star, color=:red)
    plot!(fig_3d, xlabel="x (km)", ylabel="y (km)", zlabel="z (km)",
        title="P1 vs P2 Trajectories", legend=:outertopright)

    # Time series plots for x, y, z 
    fig_xyz = plot(layout=(3, 1), size=(800, 600))

    plot!(fig_xyz[1], df.t, df.p1_x, label="P1", lw=2, color=:blue)
    plot!(fig_xyz[1], df.t, df.p2_x, label="P2", lw=2, color=:red)
    plot!(fig_xyz[1], ylabel="x (km)", title="Position vs Time")

    plot!(fig_xyz[2], df.t, df.p1_y, label="P1", lw=2, color=:blue)
    plot!(fig_xyz[2], df.t, df.p2_y, label="P2", lw=2, color=:red)
    plot!(fig_xyz[2], ylabel="y (km)")

    plot!(fig_xyz[3], df.t, df.p1_z, label="P1", lw=2, color=:blue)
    plot!(fig_xyz[3], df.t, df.p2_z, label="P2", lw=2, color=:red)
    plot!(fig_xyz[3], xlabel="Time (s)", ylabel="z (km)")

    return fig_3d, fig_xyz

end

# Plot game 1 for sanity check 
fig_3d, fig_xyz = plot_game_csv("test/game_rv_hists/p1_mixed_p2_mixed/game_1.csv")
display(fig_3d)
display(fig_xyz)

