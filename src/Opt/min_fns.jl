using LinearAlgebra 
using Optim 

## ============================================ ##

"Minimize Δv for trajectory with N segments "
function min_Δv(  
    rv_0,               # initial position vector 
    rv_f,               # final position vector 
    tof,                # time of flight 
    N      = 20,        # number of segments 
    mu     = 1.0,       # gravitational parameter 
    dm     = "pro",     # direction of motion
    Δv_max = 2.0,       # maximum Δv 
) 

    tof_N, Δv_vec = lambert_init_guess( rv_0, rv_f, tof, N, mu, dm ) 
    x_0 = reshape( Δv_vec, N*3, 1 ) 
    
    # define objective function 
    obj_fn(x) = sum_norm_Δv( x, N ) 
    
    # equality constraint 
    c_fn(x) = miss_distance_prop_kepler_Nseg( rv_0, x, N, rv_f, tof_N, mu ) 
    
    # inequality constraint ? 
    h_fn(x) = constrain_Δv( x, N, Δv_max ) 
    
    # minimize constrained 
    x_min  = min_aug_L( obj_fn, x_0, c_fn, h_fn ) 
    
    # get solution 
    Δv_sol = reshape( x_min, N, 3 ) 

    return Δv_sol 
end 

export min_Δv 

## ============================================ ##

"Minimize Δv for trajectory with N segments "
function min_Δv_dist(  
    rv_0,               # initial position vector 
    rv_f,               # final position vector 
    tof,                # time of flight 
    N      = 20,        # number of segments 
    mu     = 1.0,       # gravitational parameter 
    dm     = "pro",     # direction of motion
    Δv_max = 2.0,       # maximum Δv 
) 

    tof_N, Δv_vec = lambert_init_guess( rv_0, rv_f, tof, N, mu, dm ) 
    x_0 = reshape( Δv_vec, N*3, 1 ) 
    
    # define objective function 
    obj_fn(x) = + sum_norm_Δv( x, N ) + 
                miss_distance_prop_kepler_Nseg( rv_0, x, N, rv_f, tof_N, mu ) 
    
    # inequality constraint ? 
    h_fn(x) = constrain_Δv( x, N, Δv_max ) 
    
    # minimize constrained 
    x_min  = min_aug_L( obj_fn, x_0, nothing, h_fn ) 
    # x_min  = min_aug_L( obj_fn, x_0 ) 
    
    # get solution 
    Δv_sol = reshape( x_min, N, 3 ) 

    return Δv_sol 
end 

export min_Δv_dist 

## ============================================ ##

"Maximize Δv for trajectory with N segments (seems to be not working)"
function max_Δv_dist(  
    rv_0,               # initial position vector 
    rv_f,               # final position vector 
    tof,                # time of flight 
    N      = 20,        # number of segments 
    mu     = 1.0,       # gravitational parameter 
    dm     = "pro",     # direction of motion
    Δv_max = 2.0,       # maximum Δv 
) 

    tof_N, Δv_vec = lambert_init_guess( rv_0, rv_f, tof, N, mu, dm ) 
    x_0 = reshape( Δv_vec, N*3, 1 ) 
    
    # define objective function 
    obj_fn(x) = - sum_norm_Δv( x, N ) - 
                miss_distance_prop_kepler_Nseg( rv_0, x, N, rv_f, tof_N, mu ) 
    
    # inequality constraint ? 
    h_fn(x) = constrain_Δv( x, N, Δv_max ) 
    
    # minimize constrained 
    x_min  = min_aug_L( obj_fn, x_0, nothing, h_fn ) 
    
    # get solution 
    Δv_sol = reshape( x_min, N, 3 ) 

    return Δv_sol 
end 

export max_Δv_dist 

## ============================================ ##

"Minimize function using Optim"
function min_optim(  
    fn,                     # objective function 
    x_0,                    # initial guess 
    method = NelderMead(), 
    tol = 1e-6, 
) 

    # assign gradient fn 
    dfn = x -> ForwardDiff.gradient( fn, x ) 

    # minimize 
    result = optimize( fn, dfn, x_0, method ) 
    x_min  = result.minimizer 

    return x_min 
end

export min_optim 

## ============================================ ##

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
