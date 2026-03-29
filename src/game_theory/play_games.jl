

function prop_game_step(game, params, rng)

    # get the current player states 
    rv_E = game.p1_state[end].rv_chosen[end, :]
    rv_P = game.p2_state[end].rv_chosen[end, :]

    # return reference orbit for most recent step 
    t_ref_E, rv_ref_E, kep_ref_E = find_ref_orbit(game, params)

    # generate rv_ref_E_hist 
    t_ref_E_hist, rv_ref_E_hist = prop_rv_ref(kep_ref_E, params)

    # now move game forward one step 
    push!(game.tt, game.tt[end] + params.tt_step)
    push!(game.k_replan, game.k_replan[end] + 1)
    push!(game.rv_E, rv_E)
    push!(game.rv_P, rv_P)
    push!(game.t_ref_E, t_ref_E .+ t_ref_E_hist)
    push!(game.rv_ref_E, rv_ref_E_hist)

    # propagate chosen trajectories for evader and pursuer 
    t_E_hist, rv_E_hist, t_P_hist, rv_P_hist = prop_chosen_rv(rv_E, rv_P, params)

    # save player state and control hists (carry forward beliefs from previous step)
    p1 = player_struct([], [], [], [], [], [], rv_E_hist, [], [], [], deepcopy(game.p1_state[end].fp_belief), deepcopy(game.p1_state[end].strategy_belief), deepcopy(game.p1_state[end].tracked_strategies))
    p2 = player_struct([], [], [], [], [], [], rv_P_hist, [], [], [], deepcopy(game.p2_state[end].fp_belief), deepcopy(game.p2_state[end].strategy_belief), deepcopy(game.p2_state[end].tracked_strategies))
    players = [p1, p2]

    # compute all possible Δv solutions - 
    players = compute_states_nash(params, game, players, rng)

    # save player state in game 
    push!(game.p1_state, players[1])
    push!(game.p2_state, players[2])

    return game
end

export prop_game_step

## ====================================================================

function run_game(rng, N_replan=10, p1_strategy="mixed", p2_strategy="mixed")

    params, players, game = init_game(rng, p1_strategy, p2_strategy)

    for ii = 1:N_replan-1
        println("step: ", ii + 1, "\n")
        game = prop_game_step(game, params, rng)
    end

    return game, params
end

export run_game


## ====================================================================

function run_MC_games(rng, N_games, N_replan, p1_strategy="mixed", p2_strategy="mixed")

    games_vec = []

    for jj = 1:N_games

        params, players, game = init_game(rng, p1_strategy, p2_strategy)

        for ii = 1:N_replan-1
            println("game: ", jj, " step: ", ii + 1, "\n")
            game = prop_game_step(game, params, rng)
        end

        push!(games_vec, game)

    end

    # print_MC_stats( games_vec ) 
    save_games_vec(games_vec)

    return games_vec
end

export run_MC_games

# ==================================================================== 

using Base.Threads

function run_MC_games_parallel(rng, N_games, k_replan, p1_strategy="mixed", p2_strategy="mixed")

    # Create thread-safe storage
    games_vec = Vector{Any}(undef, N_games)

    # Create independent RNGs for each game (thread-safe)
    seeds = rand(rng, UInt, N_games)

    # run games in parallel
    Threads.@threads for i_game = 1:N_games

        local_rng = MersenneTwister(seeds[i_game])

        params, players, game = init_game(local_rng, p1_strategy, p2_strategy)

        # propogate each game forward one step at a time 
        for i_step = 1:k_replan-1
            println("game: ", i_game, " step: ", i_step + 1)
            game = prop_game_step(game, params, local_rng)
        end

        games_vec[i_game] = game

    end

    save_games_vec(games_vec)

    return games_vec
end

export run_MC_games_parallel



## ====================================================================

# using Infiltrator

# function rv_E_P_strategy(game, params)

#     # get most recent player states 
#     p1 = game.p1_state[end]
#     p2 = game.p2_state[end]

#     # p1_chosen = p1.chosen 
#     # p2_chosen = p2.chosen 

#     # choose the trajectory based on the strategy 
#     # p1_chosen, p2_chosen = p_strategy( game, game.k_replan[end], params.strategy ) 

#     # get the rv for each player at the k_tt_replan + 1 time step --> make it CURRENT state 
#     # rv_E = p1.X[ p1_chosen ][ params.k_tt_replan + 1, : ] 
#     # rv_P = p2.X[ p2_chosen ][ params.k_tt_replan + 1, : ] 

#     # @infiltrate 

#     rv_E = p1.rv_chosen[end, :]
#     rv_P = p2.rv_chosen[end, :]

#     return rv_E, rv_P
# end

# export rv_E_P_strategy

## ==================================================================== 

