## ============================================ ##

"Propagates an initial state through a vector of N trajectory segments" 
function prop_2Body_dt_Nseg(
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
        rv_k_dv = apply_dv( rv_k, Δv_vec[i,:] ) 

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
 
    return rv_hist, t_hist 
end

export prop_2Body_dt_Nseg 

#============================================================

APPLY_DV:

Description: Helper function for prop_stateUV_Nseg() and prop_stateUV_Nseg_range() that adds Δv to a state vector's velocity

Inputs:
    1. rv  - state vector
    2. Δv - velocity vector

Outputs:
    1. rv - updated state vector

============================================================#

export apply_dv 
function apply_dv(
    rv,         # state vector 
    Δv,         # velocity vector 
)

    # Applying Δv
    r = rv[1:3]
    v = rv[4:6] + Δv
    rv = vcat(r, v)

    return rv
end 


#============================================================
MISS_DISTANCE: 

Description: Calculates miss distance between trajectories

Inputs: 
    rv_0      - initial state vector of form [r; v]
    Δv_vec  - [N,3] matrix of Δv vectors, Δv_i at [i,:] 
    N       - number of segments 
    xf      - target state vector of form [r; v] 
    tof_N   - tof for each segment 
    mu      - Graviational parameter (default = 1.0) 

Outputs: 
    Δxf     - Miss distance
============================================================#

export miss_distance 
function miss_distance( 
    rv_0, 
    Δv_vec, 
    N, 
    xf, 
    tof_N = 1.0, 
    mu = 1.0, 
)

    # Propagating To Final State
    rv_hist, Δt = prop_2Body_dt_Nseg( rv_0, Δv_vec, N, tof_N, mu ) 

    # extract final state 
    xf_prop = rv_hist[end,:] 

    # Finding Miss Distance
    Δxf = abs.(xf_prop[1:3] - xf[1:3])

    return Δxf
end

