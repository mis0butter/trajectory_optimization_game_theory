using trajectory_optimization_game_theory 
using ForwardDiff 
using LinearAlgebra 

## ============================================ ##
# BFGS 

using LinearAlgebra 

# obj fn and gradient 
fn(x)  = x[1]^4 + x[1]*x[2] + 3*x[2]^2 ;
dfn(x) = [ 4*x[1]^3 + x[2] ; x[1] + 6*x[2] ]  
    
# initial guess 
x0 = [3.0, 3.0] ;

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
    xnew = x + alpha * pk ; 

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
    rhok   = 1/( yk' * sk ) ; 
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

function fn_xy( x,y )
    return fn( [x,y] )
end 

x, y = collect( -2 : 0.1 : 2 ), collect( -2 : 0.1 : 2 ) 
z = fn_xy.( x, y' ) 

x_sol  = x_hist[end] 
x_bfgs = min_bfgs( fn, dfn, x0 )

fig = plot_surface( x, y, z ) 
fig = plot_scatter( x_sol[1], x_sol[2] , fn(x_sol), fig ) 


