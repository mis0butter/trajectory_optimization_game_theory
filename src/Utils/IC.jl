using LinearAlgebra

## ====================================================================

function init_game(
    rng,
    p1_strategy="mixed",
    p2_strategy="mixed"
)

    # orbit params 
    mu = 398600.4415   # gravitational parameter 
    r = 6378.0        # Earth radius [km] 

    # radius of polygon circle!!!! 
    R_polygon = 6378.0 / 100

    # orbital elements
    a = r + 620.0
    kep0_E = [a, 0.01, 20 * pi / 180, 10.0 * pi / 180, 20.0 * pi / 180, 25.0 * pi / 180]
    rv_0_E = kep2cart(kep0_E, mu)

    # initial conditions for pursuer 
    # r_0_P = rand_IC( rv_0_E, R_polygon, rng ) 
    # rv_0_P = [ r_0_P ; rv_0_E[4:6] ]
    kep0_P = [a * 1.005, 0.01, 20 * pi / 180, 10.0 * pi / 180, 20.0 * pi / 180, 25.0 * pi / 180]
    rv_0_P = kep2cart(kep0_P, mu)

    # get period of orbit 
    T = 2 * pi * sqrt(a^3 / mu)

    # --- orbit discretization ---
    n_seg_orbit         = 50    # total segments per orbit
    n_seg_horizon       = 10    # segments per optimization window
    n_seg_per_game_step = 5     # segments executed per game step

    # --- derived time quantities ---
    dt_seg        = T / n_seg_orbit
    t_horizon     = n_seg_horizon * dt_seg
    tt_game_step  = n_seg_per_game_step * dt_seg

    # start time 
    tt = 0

    # propagate ref orbits 
    t_E, rv_E_hist = propagate_2Body(rv_0_E, t_horizon, mu, 1.0)
    t_P, rv_P_hist = propagate_2Body(rv_0_P, t_horizon, mu, 1.0)
    rv_E_hist = vv2m(rv_E_hist)
    rv_P_hist = vv2m(rv_P_hist)

    # save reference orbit 
    kep0_ref_E = copy(kep0_E)

    # game parameters 
    params = (
        mu=mu,
        r=r,
        n_seg_horizon=n_seg_horizon,
        t_horizon=t_horizon,
        n_seg_per_game_step=n_seg_per_game_step,
        tt_game_step=tt_game_step,
        kep0_ref_E=kep0_ref_E,
        T=T,
        strategy=p1_strategy,
        p2_strategy=p2_strategy,
        R_polygon=R_polygon
    )

    # save rv_ref from E to position for vertices of polygon 
    # rv_ref_E = rv_E_hist
    t_ref_E, rv_ref_E = prop_kepler_tof_Nseg(rv_0_E, zeros(params.n_seg_horizon, 3), params.n_seg_horizon, params.t_horizon / params.n_seg_horizon, params.mu)

    # save player state and control hists 
    tracked_strategies = ["mixed", "greedy", "random", 1, 2, 3, 4, 5, 6]
    strategy_belief = ones(length(tracked_strategies)) / length(tracked_strategies)
    p = player_struct([], [], [], [], [], [], [], [], [], [], ones(6), strategy_belief, tracked_strategies)
    players = [p, deepcopy(p)]

    players[1].rv_0_hist = rv_E_hist
    players[2].rv_0_hist = rv_P_hist

    # init game 
    game = game_struct([], [], [], [], [], [], [], [], [])

    push!(game.tt, tt)
    push!(game.i_game_step, 1)
    push!(game.rv_E, rv_0_E)
    push!(game.rv_P, rv_0_P)
    push!(game.t_ref_E, t_ref_E)
    push!(game.rv_ref_E, rv_ref_E)
    push!(game.params, params)

    # compute all possible Δv solutions 
    players = compute_states_nash(params, game, players, rng)

    # save player state in game 
    push!(game.p1_state, players[1])
    push!(game.p2_state, players[2])

    return params, players, game
end

export init_game


## ====================================================================

"Create dummy IC for lambert transfer and then breaking up into smaller Δv. for testing purposes only!!!"
function lambert_IC()

    # define IC, target state, and lambert solve 
    R = 6378.0                # Earth radius [km] 
    r_0 = [20.0e6, 20.0e6, 0]  # [km] 
    r_f = [-20.0e6, 10.0e6, 0]  # [km] 
    # r_0 = [ R + 500, R + 500, 0 ]
    # r_f = [ R - 500, R - 500, 0 ] 
    tof = 1.0 * 86400
    mu = 398600.4418e9         # [km^3/s^2] 
    dm = "pro"

    # solve lambert orbit 
    v_0, v_f = lambertbattin(r_0, r_f, mu, dm, tof)

    # # # let's non-dimensionalize everything 
    # r_0, v_0 = nondim_rv( r_0, v_0, mu, R )
    # r_f, v_f = nondim_rv( r_f, v_f, mu, R )  
    # mu = 1.0 

    rv_0 = [r_0; v_0]

    # propagate lambert orbit 
    t_lambert, rv_lambert = propagate_2Body(rv_0, tof, mu, 1.0)
    rv_lambert = mapreduce(permutedims, vcat, rv_lambert)

    # N segments 
    N = 20
    # Δv_vec = zeros(N, 3) 
    # Δv_vec[1,:] = v_0 

    # initial position has all z velocity 
    v_z = [0; 0; norm(v_0)]

    # compute delta v vec for lambert solution 
    dv = v_0 - v_z
    dv = v_0

    # break up delta v into smaller segments 
    Δv_vec = []
    for i = 1:N
        push!(Δv_vec, dv / N)
    end
    Δv_vec = mapreduce(permutedims, vcat, Δv_vec)

    return r_0,            # initial position 
    r_f,            # final position 
    v_0,            # initial velocity 
    v_f,            # final velocity 
    rv_lambert,     # lambert orbit 
    Δv_vec,         # delta v vector 
    tof,            # time of flight 
    N,              # number of segments 
    mu              # gravitational parameter  
end

export lambert_IC
# r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 

