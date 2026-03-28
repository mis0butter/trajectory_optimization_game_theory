function prop_opt_Nseg(
    rv_0,
    rv_f,
    tof,
    dm,
    N=10,
    mu=1.0
)

    # first compute lambert 
    _, Δv = prop_lambert_soln(rv_0, rv_f, tof, dm, mu)

    # set initial guess 
    tof_N = tof / N / 4

    Δv_vec = [Δv]
    for i = 1:N-1
        push!(Δv_vec, zeros(3))
    end
    Δv_vec = vv2m(Δv_vec)
    Δv_vec_flat = reshape(Δv_vec, N * 3, 1)
    x_0 = 0.9 * [tof_N; Δv_vec_flat]

    # define objective function 
    obj_fn(x) = sum_norm_Δv(x, N)
    obj_fn(x_0)

    # equality constraint 
    c_fn(x) = miss_distance_prop_kepler_Nseg(rv_0, x[2:end], N, rv_f, x[1], mu)
    c_fn(x_0)

    # inequality constraint ? 
    Δv_max = 2.0
    h_fn(x) = constrain_Δv(x, N, Δv_max)
    h_fn(x_0)

    # minimize constrained 
    x_min = min_aug_L(obj_fn, x_0, c_fn, h_fn)

    # get solution 
    Δv_sol = reshape(x_min[2:end], N, 3)
    tof_N_sol = x_min[1]

    return tof_N_sol, Δv_sol
end

## ====================================================================

"Propagates an initial state through a vector of N trajectory segments using dynamics integration"
function prop_2Body_tof_Nseg(
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    tof_N,          # tof for each segment 
    mu=1.0        # gravitational parameter 
)

    # Creating Iteration Variables
    rv_k = copy(rv_0)
    Δt = N * tof_N

    # Propagating Through Each Segment 
    rv_hist = [rv_k]
    t_hist = [0]
    for i = 1:N

        # apply dv 
        rv_k_dv = apply_Δv(rv_k, Δv_vec[i, :])

        # propagate and save 
        t, rv = propagate_2Body(rv_k_dv, tof_N, mu)
        for j = 2:length(rv)
            push!(rv_hist, rv[j])
        end
        t_hist = [t_hist; t_hist[end] .+ t[2:end]]

        # set up next iter 
        rv_k = rv[end]

    end
    rv_hist = mapreduce(permutedims, vcat, rv_hist)

    return t_hist, rv_hist
end

export prop_2Body_tof_Nseg


## ====================================================================

"Adds Δv to a state vector's velocity"
function apply_Δv(
    rv,         # state vector 
    Δv,         # velocity vector 
)

    # Applying Δv
    r = rv[1:3]
    v = rv[4:6] + Δv
    rv = vcat(r, v)

    return rv
end

export apply_Δv

## ====================================================================

function find_ref_orbit(game, params)

    # rv_E, rv_P = rv_E_P_strategy( game, params ) 

    t_ref_E_hist = game.t_ref_E[end]
    rv_ref_E_hist = game.rv_ref_E[end]

    # get the rv for each player at the k_tt_replan + 1 time step --> make it CURRENT state 
    t_ref_E = t_ref_E_hist[params.k_tt_replan+1, :]
    rv_ref_E = rv_ref_E_hist[params.k_tt_replan+1, :]

    kep_ref_E = cart2kep(rv_ref_E, params.mu)

    return t_ref_E, rv_ref_E, kep_ref_E
end

export find_ref_orbit


## ====================================================================

function prop_rv_ref(kep_ref_E, params)

    # save OG reference orbit 
    kep0_ref_E = params.kep0_ref_E
    kep0_ref_E[end] = kep_ref_E[end]

    rv0_ref_E = kep2cart(kep0_ref_E, params.mu)

    # t_ref_E, rv_ref_E_hist = propagate_2Body(rv0_ref_E, params.tof, params.mu, 1.0) 
    t_ref_E, rv_ref_E_hist = prop_kepler_tof_Nseg(rv0_ref_E, zeros(params.N, 3), params.N, params.tof / params.N, params.mu)
    # rv_ref_E_hist = vv2m(rv_ref_E_hist) 

    return t_ref_E, rv_ref_E_hist
end

export prop_rv_ref
