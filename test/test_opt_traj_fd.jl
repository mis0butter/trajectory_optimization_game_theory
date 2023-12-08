using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences
using LinearAlgebra 

## ============================================ ##
# test miss distance (min_bfgs exploded) 

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
rv_0  = [ r_0 ; v_0 ] 
rv_f  = [ r_f ; v_f ] 
tof_N = tof / N 

# x starts off as a [N*3, 1] vector 
Δv_vec_flat = reshape( Δv_vec, N*3, 1 ) 

function miss_Δv_flat( 
    rv_0, Δv_vec_flat, N, rv_f, tof_N, mu )
    
    Δv_vec = reshape( Δv_vec_flat, N, 3 ) 
    miss   = miss_distance_prop_kepler_Nseg( rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
    return miss 
end 

# test 
out = miss_Δv_flat( rv_0, Δv_vec_flat, N, rv_f, tof_N, mu )

# define objective fn 
fn( Δv_vec_flat ) = miss_Δv_flat( rv_0, Δv_vec_flat, N, rv_f, tof_N, mu ) 
fn( Δv_vec_flat ) 

# create gradient fn 
dfn_fd  = Δv_vec_flat -> ForwardDiff.gradient( fn, Δv_vec_flat ) 
dfn_fdm = Δv_vec_flat -> grad( central_fdm(5, 1), fn, Δv_vec_flat )[1]

# compute gradient 
g_fd  = dfn_fd( Δv_vec_flat ) 
g_fdm = dfn_fdm( Δv_vec_flat ) 
println( "err norm = ", norm( g_fd - g_fdm ) ) 

# minimize 
x_min_fd   = min_bfgs( fn, dfn_fd, Δv_vec_flat )  
x_min_fdm  = min_bfgs( fn, dfn_fdm, Δv_vec_flat )  
Δv_sol_fd  = reshape(x_min_fd, N, 3) 
Δv_sol_fdm = reshape(x_min_fdm, N, 3) 

# plot solutions 
t_kep, rv_kep_fd  = prop_kepler_tof_Nseg( rv_0, Δv_sol_fd, N, tof_N, mu ) 
t_kep, rv_kep_fdm = prop_kepler_tof_Nseg( rv_0, Δv_sol_fdm, N, tof_N, mu ) 
fig = plot_orbit( rv_kep_fd ) 
fig = plot_orbit( rv_kep_fdm, fig )  

## ============================================ ##
# test miss distance with tof_N part of obj fn 

tof_Δv = [ tof_N ; Δv_vec_flat ] 

function miss_tof_Δv_flat( 
    rv_0, 
    tof_N_Δv_vec_flat, 
    N, 
    rv_f, 
    mu 
) 

    # get tof and Δv 
    tof_N       = tof_N_Δv_vec_flat[1] 
    Δv_vec_flat = tof_N_Δv_vec_flat[2:end] 
    
    Δv_vec = reshape( Δv_vec_flat, N, 3 ) 
    miss   = miss_distance_prop_kepler_Nseg( rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
    return miss 
end 

# define objective fn 
fn( tof_N_Δv_vec_flat ) = miss_tof_Δv_flat( rv_0, tof_N_Δv_vec_flat, N, rv_f, mu ) 
fn( tof_Δv ) 

# create gradient fn  
dfn_fd  = tof_Δv -> ForwardDiff.gradient( fn, tof_Δv ) 
dfn_fdm = tof_Δv -> grad( central_fdm(5, 1), fn, tof_Δv )[1]

# compute gradient 
g_fd  = dfn_fd( tof_Δv ) 
g_fdm = dfn_fdm( tof_Δv ) 
println( "err norm = ", norm( g_fd - g_fdm ) ) 

# minimize 
x_min_fd   = min_bfgs( fn, dfn_fd, tof_Δv )  
# x_min_fdm  = min_bfgs( fn, dfn_fdm, tof_Δv )  

# get soln 
tof_N_fd  = x_min_fd[1] 
Δv_sol_fd = reshape( x_min_fd[2:end], N, 3 ) 

# plot solution 
t_kep, rv_kep_fd = prop_kepler_tof_Nseg( rv_0, Δv_sol_fd, N, tof_N, mu ) 
fig = plot_orbit( rv_kep_fd ) 

## ============================================ ##
## ============================================ ##
# now test miss distance while changing obj fn 

tof_Δv = [ tof_N ; Δv_vec_flat ] 

function miss_mag_tof_Δv_flat( 
    rv_0, 
    tof_N_Δv_vec_flat, 
    N, 
    rv_f, 
    mu 
) 

    # get tof and Δv 
    tof_N       = tof_N_Δv_vec_flat[1] 
    Δv_vec_flat = tof_N_Δv_vec_flat[2:end] 
    
    Δv_vec = reshape( Δv_vec_flat, N, 3 ) 
    miss   = miss_distance_prop_kepler_Nseg( rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

    # state magnitude 
    state_mag = norm( tof_N_Δv_vec_flat ) 
    
    return miss + state_mag 
end 

# define objective fn 
fn( tof_N_Δv_vec_flat ) = miss_mag_tof_Δv_flat( rv_0, tof_N_Δv_vec_flat, N, rv_f, mu ) 
fn( tof_Δv ) 

# create gradient fn  
dfn_fd  = tof_Δv -> ForwardDiff.gradient( fn, tof_Δv ) 
# dfn_fdm = tof_Δv -> grad( central_fdm(5, 1), fn, tof_Δv )[1] 

# compute gradient 
g_fd  = dfn_fd( tof_Δv ) 
# g_fdm = dfn_fdm( tof_Δv ) 
# println( "err norm = ", norm( g_fd - g_fdm ) ) 

# minimize 
x_min_fd   = min_bfgs( fn, dfn_fd, tof_Δv )  
# x_min_fdm  = min_bfgs( fn, dfn_fdm, tof_Δv )  

# get soln 
tof_N_fd  = x_min_fd[1] 
Δv_sol_fd = reshape( x_min_fd[2:end], N, 3 ) 

# plot solution 
t_kep, rv_kep_fd = prop_kepler_tof_Nseg( rv_0, Δv_sol_fd, N, tof_N_fd, mu ) 
fig = plot_orbit( rv_kep_fd ) 


 


