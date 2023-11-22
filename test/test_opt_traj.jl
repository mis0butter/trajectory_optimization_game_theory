using trajectory_optimization_game_theory 
using ForwardDiff 
using LinearAlgebra 

## ============================================ ##
# test miss distance 

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
rv_0 = [ r_0 ; v_0 ] 
rv_f = [ r_f ; v_f ] 
tof_N = tof / N 

miss_kepler = miss_distance_prop_kepler( 
    rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

obj( Δv_vec ) = miss_distance_prop_kepler( 
        rv_0, Δv_vec, N, rv_f, tof_N, mu ) 

dobj = Δv_vec -> ForwardDiff.jacobian( obj, Δv_vec ) 

Δv_soln = min_bfgs( obj, dobj, Δv_vec )

# using Optim 

# # optimization 
# od       = OnceDifferentiable( obj, Δv_vec ; autodiff = :forward ) 
# result   = optimize( od, Δv_vec, BFGS() ) 
# hp       = result.minimizer 


## ============================================ ##
# test miss distance (min_bfgs exploded) 

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
rv_0 = [ r_0 ; v_0 ] 
rv_f = [ r_f ; v_f ] 
tof_N = tof / N 

# x starts off as a [N*3, 1] vector 
Δv_vec_flat = reshape( Δv_vec, N*3, 1 ) 

function miss_Δv_flat( rv_0, Δv_vec_flat, N, rv_f, tof_N, mu )
    
    Δv_vec = reshape( Δv_vec_flat, N, 3 ) 
    miss   = miss_distance_prop_kepler( 
    rv_0, Δv_vec, N, rv_f, tof_N, mu ) 
    
    return miss 
end 

out = miss_Δv_flat( rv_0, Δv_vec_flat, N, rv_f, tof_N, mu )

obj_fn( Δv_vec_flat ) = miss_Δv_flat( rv_0, Δv_vec_flat, N, rv_f, tof_N, mu )

obj_fn( Δv_vec_flat ) 
dobj_fn = Δv_vec_flat -> ForwardDiff.gradient( obj_fn, Δv_vec_flat ) 

dobj_fn( Δv_vec_flat ) 

## ============================================ ##
    
    
# initial guess 
# x0 = [3.0, 3.0] ;
x0 = Δv_vec 

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
    while fn(xnew) > fn(x) + c * alpha * g' * pk 
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
    # rhok   = 1/( yk' * sk ) ; 
    rhok   = inv( yk' * sk ) 
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
