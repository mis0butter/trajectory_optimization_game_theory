using LinearAlgebra 

"Minimize a function using BFGS method"
function min_bfgs(  
    fn,                 # objective function 
    dfn,                # gradient of objective function 
    x0,                 # initial guess 
    tol     = 1e-6,     # termination tolerance 
    maxiter = 1000,     # maximum number of allowed iterations 
    dxmin   = 1e-6,     # minimum allowed perturbation 
    beta    = 0.707,    # backtracking line search parameter 
    c       = 1e-4,     # backtracking line search parameter 
) 

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

        # search direction 
        g    = dfn(x)               # gradient   
        pk   = - Hk * g             # search direction 
        xnew = x + alpha .* pk      # take step 
    
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
    
        # secant equation - update Hk  
        Hk_new = update_Hk( x, xnew, Hk, dfn ) 
    
        # update termination metrics
        dx    = norm(xnew-x) ; 
        x     = xnew ;
        Hk    = Hk_new ; 
        niter = niter + 1 ;
        if niter == maxiter 
            println("maxiter exceeded") 
        end
    
        # save hist 
        push!( x_hist, x )
        push!( f_hist, fn(x) )  
    
    end 

    x_sol = x_hist[end] 

    return x_sol 
end 

export min_bfgs 

## ============================================ ##

function min_golden_ratio(
    x, 
    fn; 
    tol = 1e-8, 
    itermax = 50, 
)

    # Initializing
    ϕ  = (3 - √5)/2
    xL = x[1]            ; fL = f(xL)
    xR = x[3]            ; fR = f(xR)
    x1 = xL + ϕ*(xR - xL); f1 = f(x1)
    x2 = xR - ϕ*(xR - xL); f2 = f(x2)
    δx = Inf
    iter = 0

    # Iterating
    while abs(δx) > tol
        iter += 1

        # Branching
        if f1 > f2
            # Adjusting Bounded Points
            xL = x1; fL = f1;
            x1 = x2; f1 = f2;
            
            # Re-Evaluating New Points
            x2 = xR - ϕ*(xR - xL);
            f2 = f(x2);
        else
            # Adjusting Bounded Points
            xR = x2; fR = f2;
            x2 = x1; f2 = f1;
            
            # Re-Evaluating New Points
            x1 = xL + ϕ*(xR - xL);
            f1 = f(x1);
        end

        # Finding Bounds
        δx = norm(xR - xL)

        if iter == itermax; 
            @warn "GoldRatioMin iteration failure"
            break
        end
    end

    # Finding Best Point
    X = [xL, x1, x2, xR]
    Z = [fL, f1, f2, fR]
    idx = findmin(Z)[2]

    xf = X[idx]

    return xf, iter
end


## ============================================ ##

"Update approximate Hessian Hk via BFGS secant equation (Nocedal)"
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

export update_Hk 

## ============================================ ##

"Update approximate Hessian Hk via BFGS secant equation (Russell)"
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

export update_Qk 
