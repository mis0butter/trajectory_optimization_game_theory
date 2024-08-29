using LinearAlgebra 

## ============================================ ##

function init_game( rng, strategy = "mixed" ) 

    # orbit params 
    mu  = 398600.4415   # gravitational parameter 
    r   = 6378.0        # Earth radius [km] 
    
    # orbital elements
    a      = r + 620.0  
    kep0_E = [ a, 0.01, 20*pi/180, 10.0*pi/180, 20.0*pi/180, 25.0*pi/180 ]
    rv_0_E = kep2cart(kep0_E, mu) 

    # initial conditions for pursuer 
    r_0_P = rand_IC( rv_0_E ) 
    rv_0_P = [ r_0_P ; rv_0_E[4:6] ]

    # get period of orbit 
    T = 2*pi*sqrt( a^3 / mu ) 

    # orbit parameters of pursuer and evader 
    k_tt_replan = 5       # replan every 5 * tof/N (100) seconds!!! 
    tof = T / k_tt_replan # tof for pursuer to catch up to evader  
    N   = 10              # segments 
    
    # start time 
    tt = 0 
    tt_step = k_tt_replan * tof / N 
    
    # propagate ref orbits 
    t_E, rv_E_hist = propagate_2Body(rv_0_E, tof, mu, 1.0) 
    t_P, rv_P_hist = propagate_2Body(rv_0_P, tof, mu, 1.0) 
    rv_E_hist = vv2m(rv_E_hist) 
    rv_P_hist = vv2m(rv_P_hist) 

    # save reference orbit 
    kep0_ref_E = copy( kep0_E ) 
    
    # game parameters 
    params = ( mu = mu, r = r, N = N, tof = tof, k_tt_replan = k_tt_replan, tt_step = tt_step, kep0_ref_E = kep0_ref_E, T = T, strategy = strategy ) 
    
    # save rv_ref from E to position for vertices of polygon 
    rv_ref_E_polygon = rv_E_hist 
    
    # save player state and control hists 
    p = player_struct( [], [], [], [], [], [], [] ) 
    players = [ p, deepcopy(p) ]  
    
    players[1].rv_0_hist = rv_E_hist 
    players[2].rv_0_hist = rv_P_hist 
    
    # init game 
    game = game_struct( [], [], [], [], [], [], [] ) 
    
    push!( game.tt, tt ) 
    push!( game.k_replan, 1 )
    push!( game.rv_E, rv_0_E ) 
    push!( game.rv_P, rv_0_P ) 
    push!( game.rv_ref_E_polygon, rv_ref_E_polygon )  

    # compute all possible Δv solutions 
    players = players_states( params, game, players, rng ) 

    # save player state in game 
    push!( game.p1_state, players[1] ) 
    push!( game.p2_state, players[2] )  


    return params, players, game 
end 

export init_game 


## ============================================ ##

"Create dummy IC for lambert transfer and then breaking up into smaller Δv. for testing purposes only!!!"
function lambert_IC() 

    # define IC, target state, and lambert solve 
    R   = 6378.0                # Earth radius [km] 
    r_0 = [ 20.0e6, 20.0e6, 0]  # [km] 
    r_f = [-20.0e6, 10.0e6, 0]  # [km] 
    # r_0 = [ R + 500, R + 500, 0 ]
    # r_f = [ R - 500, R - 500, 0 ] 
    tof = 1.0 * 86400 
    mu  = 398600.4418e9         # [km^3/s^2] 
    dm  = "pro" 

    # solve lambert orbit 
    v_0, v_f = lambertbattin( r_0, r_f, mu, dm, tof ) 

    # # # let's non-dimensionalize everything 
    # r_0, v_0 = nondim_rv( r_0, v_0, mu, R )
    # r_f, v_f = nondim_rv( r_f, v_f, mu, R )  
    # mu = 1.0 

    rv_0     = [ r_0; v_0 ] 

    # propagate lambert orbit 
    t_lambert, rv_lambert = propagate_2Body( rv_0, tof, mu, 1.0 ) 
    rv_lambert = mapreduce( permutedims, vcat, rv_lambert ) 

    # N segments 
    N = 20 
    # Δv_vec = zeros(N, 3) 
    # Δv_vec[1,:] = v_0 

    # initial position has all z velocity 
    v_z = [ 0; 0; norm(v_0) ]

    # compute delta v vec for lambert solution 
    dv = v_0 - v_z   
    dv = v_0 

    # break up delta v into smaller segments 
    Δv_vec = [] 
    for i = 1 : N 
        push!( Δv_vec, dv/N ) 
    end 
    Δv_vec = mapreduce( permutedims, vcat, Δv_vec ) 

    return  r_0,            # initial position 
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

