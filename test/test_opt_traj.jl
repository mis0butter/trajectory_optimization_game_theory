using trajectory_optimization_game_theory 
using ForwardDiff 
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
    miss   = miss_distance_prop_kepler( 
        rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
    return miss 
end 

out = miss_Δv_flat( 
    rv_0, Δv_vec_flat, N, rv_f, tof_N, mu )

# define objective fn 
fn( Δv_vec_flat ) = miss_Δv_flat( 
    rv_0, Δv_vec_flat, N, rv_f, tof_N, mu ) 
fn( Δv_vec_flat ) 

# create gradient 
dfn = Δv_vec_flat -> ForwardDiff.gradient( 
    fn, Δv_vec_flat ) 
dfn( Δv_vec_flat ) 

## ============================================ ##
# BFGS 
    
# initial guess 
# x0 = [3.0, 3.0] ;
x0 = Δv_vec_flat 

# loop options 
tol     = 1e-6 ;   # termination tolerance
maxiter = 1000 ;   # maximum number of allowed iterations
dxmin   = 1e-6 ;   # minimum allowed perturbation

# extract + define step size 
beta  = 0.707 ; c = 1e-4 ; 
alpha = 1     ; alpha0 = alpha ; 

# initialize gradient norm, optimization vector, iteration counter, perturbation
g = Inf ; x = x0 ; niter = 0 ; dx = Inf ;

# define secant equation (BFGS method) 
#     Bk0 = d2fn(x0) ;       Bk = Bk0 ; 
# Hk0 = inv(Bk0) ;    
Hk0 = I ; Hk = Hk0 ;  

# init hists 
x_hist = [ x0 ] ; 
f_hist = [ fn(x0) ] ; 

# BFGS 
k = 0 
while norm(g) >= tol && niter <= maxiter && dx >= dxmin 

    # increase iter 
    k += 1 
    println( "k = ", k ) 

    # calculate gradient 
    g = dfn(x) ; 

    # search direction 
    pk = - Hk * g ; 

    # take step:
    xnew = x + alpha .* pk ; 

    # backtracking line search 
    alpha = alpha0 ; 
    while fn(xnew) > fn(x) + ( c * alpha * g' * pk )[1] || isnan(fn(xnew))
        alpha = alpha * beta ; 
        xnew  = x + alpha * pk ; 
    end 

    # check step 
    if ~isfinite(norm(xnew)) 
        println("x is inf or NaN") 
    end 

    # secant equation 
    sk     = xnew - x ; 
    yk     = dfn(xnew) - dfn(x) ; 
    rhok   = 1/( yk' * sk )[1] ; 
    # rhok   = inv( yk' * sk ) 
    # I      = eye(size(Hk)) ; 
    Hk_new = ( I - rhok*sk*yk' ) * Hk * ( I - rhok*yk*sk' ) + rhok*sk*sk' ; 

    # update termination metrics
    niter = niter + 1 ;
    dx    = norm(xnew-x) ; 
    x     = xnew ;
    Hk    = Hk_new ; 

    # save function output 
    f = fn(x) ; 

    # save hist 
    push!( x_hist, x )
    push!( f_hist, f )  

end 

## ============================================ ##
# get solution 

# tof_N = tof / N / 2.4 

# define objective fn 
fn( Δv_vec_flat ) = miss_Δv_flat( 
    rv_0, Δv_vec_flat, N, rv_f, tof_N, mu ) 
fn( Δv_vec_flat ) 

# create gradient 
dfn = Δv_vec_flat -> ForwardDiff.gradient( 
    fn, Δv_vec_flat ) 
dfn( Δv_vec_flat ) 

x_min = min_bfgs( fn, dfn, Δv_vec_flat ) 

Δv_sol = reshape(x_min, N, 3) 
# Δv_sol = reshape(x_hist[end], N, 3) 

miss_kepler = miss_distance_prop_kepler( 
    rv_0, Δv_sol, N, rv_f, tof_N, mu)

t_kep, rv_kep = prop_kepler_tof_Nseg( 
    rv_0, Δv_sol, N, tof_N, mu ) 

fig = plot_orbit( rv_kep ) 

## ============================================ ##
# test tof gradient 

function miss_Δv_flat( 
    rv_0, 
    Δv_vec_flat, 
    N, 
    rv_f, 
    tof_N, 
    mu 
)
    
    Δv_vec = reshape( Δv_vec_flat, N, 3 ) 
    miss   = miss_distance_prop_kepler( 
        rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
    return miss 
end 

# define objective fn 
fn( Δv_vec_flat ) = miss_Δv_flat( 
    rv_0, Δv_vec_flat, N, rv_f, tof_N, mu ) 
fn( Δv_vec_flat ) 

# create gradient 
dfn = Δv_vec_flat -> ForwardDiff.gradient( 
    fn, Δv_vec_flat  ) 
dfn( Δv_vec_flat ) 

# define objective fn 
fn( tof_N ) = miss_Δv_flat( 
    rv_0, Δv_vec_flat, N, rv_f, tof_N, mu ) 
fn( tof_N ) 

# create gradient 
dfn = tof_N -> ForwardDiff.gradient( 
    fn, Δv_vec_flat  ) 
dfn( tof_N ) 

## ============================================ ##
# test miss distance with minimizing tof 

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
    miss   = miss_distance_prop_kepler( 
        rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
    return miss 
end 

# define objective fn 
fn( tof_N_Δv_vec_flat ) = miss_tof_Δv_flat( 
    rv_0, tof_N_Δv_vec_flat, N, rv_f, mu ) 
fn( tof_Δv ) 

# create gradient 
dfn = tof_Δv -> ForwardDiff.gradient( 
    fn, tof_Δv ) 
    
dfn( tof_Δv ) 

## ============================================ ##
# test forwarddiff gradient with finite differencing 

# define objective fn 
fn( Δv_vec_flat ) = miss_Δv_flat( 
    rv_0, Δv_vec_flat, N, rv_f, tof_N, mu ) 
fn( Δv_vec_flat ) 

# create gradient 
dfn = Δv_vec_flat -> ForwardDiff.gradient( 
    fn, Δv_vec_flat  ) 
g = dfn( Δv_vec_flat ) 

# try using finite differencing from library 
central_fdm(5, 1)(sin, 1.0) 
x = 1.0 
central_fdm(5, 1)(sin, x) 
x = [1.0, 2.0]
central_fdm(5, 1)(sin, x[1] + x[2]) 
central_fdm(5, 1)(sin, x[1] + x[2]) 

# gradient 
a    = randn(3, 3); a = a * a'
f(x) = 0.5 * x' * a * x
x    = rand(3) 
g_fdm = grad( central_fdm(5, 1), f, x )[1] 

# compare with ForwardDiff 
g_fd  = ForwardDiff.gradient( f, x )  

# test on our problem 
g_fdm = grad( central_fdm(5, 1), fn, Δv_vec_flat )[1]  
g_fd  = ForwardDiff.gradient( fn, Δv_vec_flat ) 

# print out difference 
println( "err norm: g_fdm - g_fd = ", norm(g_fdm - g_fd) ) 

dfn_fdm = Δv_vec_flat -> grad( central_fdm(5, 1), fn, Δv_vec_flat )[1] 
g_fdm = dfn_fdm( Δv_vec_flat ) 
g_fd  = dfn(Δv_vec_flat) 

# print out difference 
println( "err norm: g_fdm - g_fd = ", norm(g_fdm - g_fd) ) 

# try with BFGS 

x_min_fd  = min_bfgs( fn, dfn, Δv_vec_flat )  
x_min_fdm = min_bfgs( fn, dfn_fdm, Δv_vec_flat )  


Δv_sol_fd  = reshape(x_min_fd, N, 3) 
Δv_sol_fdm = reshape(x_min_fdm, N, 3) 

# Δv_sol = reshape(x_hist[end], N, 3) 

miss_kepler_fd  = miss_distance_prop_kepler( 
    rv_0, Δv_sol_fd, N, rv_f, tof_N, mu)
miss_kepler_fdm = miss_distance_prop_kepler( 
    rv_0, Δv_sol_fdm, N, rv_f, tof_N, mu)
    
t_kep, rv_kep_fd  = prop_kepler_tof_Nseg( 
    rv_0, Δv_sol_fd, N, tof_N, mu ) 
t_kep, rv_kep_fdm = prop_kepler_tof_Nseg( 
    rv_0, Δv_sol_fdm, N, tof_N, mu ) 
    
fig = plot_orbit( rv_kep_fd ) 
fig = plot_orbit( rv_kep_fdm, fig ) 










