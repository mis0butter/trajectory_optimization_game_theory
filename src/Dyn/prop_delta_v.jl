## ============================================ ##

"Propagates an initial state through a vector of N trajectory segments using dynamics integration" 
function prop_2Body_tof_Nseg(
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    tof_N,          # tof for each segment 
    mu = 1.0        # gravitational parameter 
) 
    
    # Creating Iteration Variables
    rv_k = copy(rv_0)
    Δt   = N * tof_N 

    # Propagating Through Each Segment 
    rv_hist = [ rv_k ] 
    t_hist  = [ 0 ] 
    for i = 1 : N 

        # apply dv 
        rv_k_dv = apply_Δv( rv_k, Δv_vec[i,:] ) 

        # propagate and save 
        t, rv = propagate_2Body( rv_k_dv, tof_N, mu ) 
        for j = 2 : length(rv) 
            push!( rv_hist, rv[j] ) 
        end 
        t_hist = [ t_hist ; t_hist[end] .+ t[2:end] ]

        # set up next iter 
        rv_k = rv[end]

    end 
    rv_hist = mapreduce( permutedims, vcat, rv_hist ) 
 
    return t_hist, rv_hist 
end

export prop_2Body_tof_Nseg 

## ============================================ ## 

"Propagate an initial state through a vector of N trajectory segments using Kepler's equations" 
function prop_kepler_tof_Nseg(
    rv_0,           # initial state vector of form [r; v] 
    Δv_vec,         # [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N,              # number of segments 
    tof_N,          # tof for each segment 
    mu = 1.0        # gravitational parameter 
) 
    
    # set up time and state hists 
    rv_hist = [ apply_Δv( rv_0, Δv_vec[1,:] ) ] 
    t_hist  = [ ] ; push!( t_hist, 0.0 ) 

    # propagate Through Each Segment 
    rv_k = copy(rv_0)
    for i = 1 : N 

        # apply dv and propagate 
        rv_k_dv = apply_Δv( rv_k, Δv_vec[i,:] ) 
        rv_k    = prop_kepler_tof( rv_k_dv, tof_N, mu ) 

        # save 
        push!( rv_hist, rv_k ) 
        push!( t_hist, t_hist[end] .+ tof_N ) 

    end 
    rv_hist = mapreduce( permutedims, vcat, rv_hist ) 
 
    return  t_hist, 
            rv_hist 
end

export prop_kepler_tof_Nseg 
# t_kep, rv_kep = prop_kepler_tof_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

## ============================================ ##

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

