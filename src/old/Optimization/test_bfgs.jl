using LinearAlgebra
using ForwardDiff

## ============================================ ##
# set-up 

# initial guess 
x0 = [ 2.0, 2.0, 2.0 ] 

# obj fn - true min at [-1, 0, 0]
fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 
dfn = x -> ForwardDiff.gradient( fn, x ) 


tol     = 1e-6     # termination tolerance 
maxiter = 1000     # maximum number of allowed iterations 
dxmin   = 1e-6     # minimum allowed perturbation 
beta    = 0.707    # backtracking line search parameter 
c       = 1e-4     # backtracking line search parameter 

## ============================================ ##

function update_Hk( x_k, x_kp1, H_k, dfn ) 

    # secant equation 
    s_k     = x_kp1 - x_k  
    y_k     = dfn(x_kp1) - dfn(x_k)   
    rho_k   = 1/( y_k' * s_k )[1]  
    # rhok   = inv( yk' * sk ) 
    # I      = eye(size(Hk)) ; 
    H_kp1 = ( I - rho_k*s_k*y_k' ) * H_k * ( I - rho_k*y_k*s_k' ) + rho_k*s_k*s_k'  

    return H_kp1 
end 

function update_Qk( x_k, x_kp1, Q_k, dfn ) 

    y = dfn(x_kp1) - dfn(x_k) 
    p = x_kp1 - x_k 
    A = Q_k * y * p' 
    τ = y' * Q_k * y 
    σ = p' * y 
    
    ΔQ_k = ( σ + τ ) / σ^2 * p * p' - 1/σ * (A + A')  

    Q_kp1 = Q_k + ΔQ_k 

    return Q_kp1 
end 

## ============================================ ##
# BFGS 

# init step size 
alpha = 1 ;     alpha0 = copy(alpha) ; 

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
while ( norm(g) >= tol ) && 
        ( niter <= maxiter ) && 
        ( dx >= dxmin ) 

    # increase iter 
    k += 1 

    # calculate gradient 
    g = dfn(x) ; 

    # search direction 
    pk = - Hk * g ; 

    # take step:
    xnew = x + alpha .* pk ; 

    # backtracking line search 
    alpha = copy(alpha0) ; 
    while fn(xnew) > fn(x) + ( c * alpha * g' * pk )[1] || 
            isnan(fn(xnew))
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

    Hk_new_test = update_Hk( x, xnew, Hk, dfn ) 
    Qk_new_test = update_Qk( x, xnew, Hk, dfn ) 

    # update termination metrics
    dx    = norm(xnew-x) ; 
    x     = xnew ;
    Hk    = Hk_new ; 
    niter = niter + 1 ;
    if niter == maxiter 
        println("maxiter exceeded") 
    end

    # save function output 
    f = fn(x) ; 

    # save hist 
    push!( x_hist, x )
    push!( f_hist, f )  

end 

x_sol = x_hist[end] 

