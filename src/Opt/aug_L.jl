using LinearAlgebra 

## ============================================ ##

"Construct Augmented Lagrangian function 
= obj_fn(x) + λ' * ψ_fn(x) + (p./2)' * ψ_fn(x).^2 "
function aug_L_fn( 
    obj_fn,     # objective function 
    ψ_fn,       # constraint function(s) 
    x,          # state vector 
    λ,          # Lagrange multiplier(s) 
    p,          # penalty parameter(s) 
) 

    augL = obj_fn(x) 
    ψ    = ψ_fn(x) 

    for i in eachindex(λ) 
        augL += λ[i] * ψ[i] + (p[i]/2) * ψ[i]^2 
    end 

    return augL 
end 

export aug_L_fn 

## ============================================ ##

"Construct inequality-constrained Augmented Lagrangian function 
= obj_fn(x) - λ^2 / 2p OR 
= obj_fn(x) + λ' * ψ_fn(x) + (p./2)' * ψ_fn(x).^2 "
function aug_L_ineq_fn( 
    obj_fn,     # objective function 
    h_fn,       # inequality constraint function(s) 
    x,          # state vector 
    λ,          # Lagrange multiplier(s) 
    p,          # penalty parameter(s) 
) 

    augL = obj_fn(x) 
    h    = h_fn(x)

    for i in eachindex(λ) 
        if h[i] < -λ[i]/p[i] 
            augL -= λ[i]^2 / (2*p[i]) 
        else
            augL += λ[i] * h[i] + (p[i]/2) * h[i]^2 
        end 
    end 

    return augL 
end 

export aug_L_ineq_fn 



## ============================================ ##

"Construct equality and inequality-constrained Augmented Lagrangian function 
= obj_fn(x) - λ^2 / 2p OR 
= obj_fn(x) + λ' * ψ_fn(x) + (p./2)' * ψ_fn(x).^2 "
function aug_L_eq_ineq_fn( 
    obj_fn,     # objective function 
    c_fn,       # equality constraint function(s) 
    h_fn,       # inequality constraint function(s) 
    x,          # state vector 
    λ_c,        # Lagrange multiplier(s) 
    λ_h,        # Lagrange multiplier(s) 
    p_c,        # penalty parameter(s) 
    p_h,        # penalty parameter(s) 
) 

    augL = obj_fn(x) 
    h    = h_fn(x)
    c    = c_fn(x) 

    # equality constraints 
    for i in eachindex(λ_c) 
        augL += λ_c[i] * c[i] + (p_c[i]/2) * c[i]^2 
    end 

    # inequality constraints 
    for i in eachindex(λ_h) 
        if h[i] < -λ_h[i]/p_h[i] 
            augL -= λ_h[i]^2 / (2*p_h[i]) 
        else
            augL += λ_h[i] * h[i] + (p_h[i]/2) * h[i]^2 
        end 
    end 

    return augL 
end 

export aug_L_eq_ineq_fn 

## ============================================ ##

"Minimize equality-constrained Augmented Lagrangian"
function min_aug_L_eq( 
    obj_fn,                                         # objective function 
    c_fn,                                           # constraint function: c = 0 
    x_0,                                            # initial guess 
    tol     = 1e-6, 
    λ_0     = zeros(length(c_fn(x_0))),             # initial Lagrange multiplier  
    p_0     = 10.0 * ones(length(c_fn(x_0))),       # initial penalty parameter  
    γ       = 2.0,                                  # penalty parameter 
) 

    # step 0: initialize 
    λ_k = copy( λ_0 )
    p_k = copy( p_0 ) 
    x_k = copy( x_0 ) 

    k = 0 ; loop = true 
    while loop 

        k += 1 

        # step 1: assign augmented Lagrangian fn 
        fn(x_k) = aug_L_fn( obj_fn, c_fn, x_k, λ_k, p_k ) 
        dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

        # step 2: minimize unconstrained problem  
        x_min = min_bfgs( fn, dfn, x_k )  
        
        # step 3 check convergence ... 
        dx = norm( x_min - x_k ) 
        if  ( dx < tol ) && 
            ( norm( c_fn(x_min) ) < tol ) 
                loop = false 
        else 
            x_k = x_min 
        end 

        # step 3: check constraint function and update parameters 
        λ_k += p_k .* c_fn(x_k) 
        p_k *= γ 

    end 

    return x_k 
end 

export min_aug_L_eq 

## ============================================ ##

"Minimize inequality-constrained Augmented Lagrangian" 
function min_aug_L_ineq( 
    obj_fn,                                         # objective function 
    h_fn,                                           # inequality constraint function: h <= 0  
    x_0,                                            # initial guess 
    tol     = 1e-6, 
    λ_0     = zeros(length(h_fn(x_0))),             # initial Lagrange multiplier  
    p_0     = 10.0 * ones(length(h_fn(x_0))),       # initial penalty parameter  
    γ       = 2.0,                                  # penalty parameter 
) 

    # step 0: initialize 
    λ_k = copy( λ_0 ) ;     p_k = copy( p_0 ) 
    x_k = copy( x_0 ) ;     h_k = h_fn( x_k )

    k = 0 ; loop = true 
    while loop 

        k += 1 

        # step 1: assign augmented Lagrangian fn 
        fn(x_k) = aug_L_ineq_fn( obj_fn, h_fn, x_k, λ_k, p_k ) 
        dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

        # step 2: minimize unconstrained problem  
        x_min = min_bfgs( fn, dfn, x_k )  

        # step 3 check convergence ... 
        dx = norm(x_min - x_k) 
        if dx < tol 
            loop = false 
        else 
            x_k = x_min 
        end 

        # update constraint values 
        h_k = h_fn( x_k )
        λ_k, p_k = update_λ_p_ineq( λ_k, p_k, h_k, γ ) 

    end 

    return x_k 
end 

export min_aug_L_ineq 

## ============================================ ##

"Minimize equality and inequality-constrained Augmented Lagrangian"
function min_aug_L_eq_ineq(  
    obj_fn,                                                         # objective function 
    c_fn,                                                           # constraint function: c = 0 
    h_fn,                                                           # inequality constraint function: h <= 0  
    x_0,                                                            # initial guess 
    tol     = 1e-6, 
    λ_0     = zeros(length( [ c_fn(x_0) ; h_fn(x_0) ] )),           # initial Lagrange multiplier  
    p_0     = 10.0 * ones(length( [ c_fn(x_0) ; h_fn(x_0) ] )),     # initial penalty parameter  
    γ       = 2.0,                                                  # penalty parameter 
) 

    # step 0: initialize 
    λ_k = copy( λ_0 ) ;     p_k = copy( p_0 ) 
    x_k = copy( x_0 ) ;     h_k = h_fn( x_k )

    # first split λ and p into eq and ineq constraints 
    N_c = length( c_fn(x_k) ) 
    λ_c = λ_k[1:N_c] ;      p_c = p_k[1:N_c] 
    λ_h = λ_k[N_c+1:end] ;  p_h = p_k[N_c+1:end] 

    k = 0 ; loop = true 
    while loop 

        # step 1: assign augmented Lagrangian fn 
        fn(x_k) = aug_L_fn( obj_fn, c_fn, x_k, λ_c, p_c ) + 
                  aug_L_ineq_fn( obj_fn, h_fn, x_k, λ_h, p_h ) 
        # fn(x_k) = aug_L_eq_ineq_fn( obj_fn, c_fn, h_fn, x_k, λ_c, λ_h, p_c, p_h )
        dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

        # step 2: minimize unconstrained problem  
        x_min = min_bfgs( fn, dfn, x_k )  

        # step 3 check convergence ... 
        dx = norm(x_min - x_k) 
        if  ( dx < tol ) && 
            ( norm(c_fn(x_min)) < tol )  
        #    ( norm(h_fn(x_min)) < tol ) 
                loop = false 
        else 
            x_k = x_min 
        end 

        # update equality constraints first 
        if norm( c_fn(x_k) ) > tol 
            λ_c += p_c .* c_fn(x_k) 
            p_c *= γ 
        end 

        # now inequality constraints 
        h_k      = h_fn(x_k) 
        λ_h, p_h = update_λ_p_ineq( λ_h, p_h, h_k, γ ) 

    end 

    return x_k
end 

export min_aug_L_eq_ineq 

## ============================================ ##

"Update Lagrange multipliers and penalty parameters for inequality constraints"
function update_λ_p_ineq( 
    λ_k,        # Lagrange multiplier(s) 
    p_k,        # penalty parameter(s) 
    h_k,        # evaluated inequality constraint function: h <= 0 
    γ           # step size 
) 

    # update parameters
    if length(λ_k) == 1 

        λ_k = max( λ_k + p_k * h_k , 0.0 ) 
        if h_k > 0 
            p_k *= γ 
        end 

    else 

        for i in eachindex(λ_k)
            λ_k[i] = max( λ_k[i] + p_k[i] * h_k[i] , 0.0 ) 
            if h_k[i] > 0 
                p_k[i] *= γ 
            end     
        end     

    end 

    return λ_k, p_k 
end 

export update_λ_p_ineq 


