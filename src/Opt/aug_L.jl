using LinearAlgebra 

## ============================================ ##

"Construct Augmented Lagrangian function 
= obj_fn(x) + λ' * ψ_fn(x) + (p./2)' * ψ_fn(x).^2 "
function aug_L_fn( 
    obj_fn,     # objective function 
    ψ_fn,       # constraint function 
    x,          # state vector 
    λ,          # Lagrange multiplier 
    p,          # penalty parameter 
) 

    augL = obj_fn(x) 
    for i in eachindex(λ)
        augL += λ[i] * ψ_fn(x)[i] + (p[i]/2) * ψ_fn(x)[i]^2 
    end 

    return augL 
end 

export aug_L_fn 

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
        dx = norm(x_min - x_k) 
        if dx < tol 
            loop = false 
        else 
            x_k = x_min 
        end 

        # step 3: check constraint function and update parameters 
        if norm( c_fn(x_k) ) > tol 
            λ_k += p_k .* c_fn(x_k) 
            p_k *= γ 
        else 
            loop = false 
        end 

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
    λ_k = copy( λ_0 )
    p_k = copy( p_0 ) 
    x_k = copy( x_0 ) 
    h_k = h_fn( x_k )

    k = 0 ; loop = true 
    while loop 
    
        k += 1 
    
        # step 1: assign augmented Lagrangian fn 
        fn(x_k) = aug_L_fn( obj_fn, h_fn, x_k, λ_k, p_k ) 
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
    
        # ... and update parameters 
        for i in eachindex(λ_k)
    
            # update λ
            λ_k[i] = max( λ_k[i] + p_k[i] * h_k[i] , 0.0 ) 
    
            # update p 
            if h_k[i] > 0 
                p_k[i] *= γ 
            end 
    
        end 
    
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
    λ_k = copy( λ_0 )
    p_k = copy( p_0 ) 
    x_k = copy( x_0 ) 
    h_k = h_fn( x_k ) 

    # get number of equality constraints 
    N_c = length( c_fn(x_0) ) 
    
    # put constraints together 
    ψ_fn(x) = [ c_fn(x) ; h_fn(x) ]  

    k = 0 ; loop = true 
    while loop 

        # assign augmented Lagrangian fn 
        fn(x_k) = aug_L_fn( obj_fn, ψ_fn, x_k, λ_k, p_k ) 
        dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

        # minimize unconstrained problem  
        x_min = min_bfgs( fn, dfn, x_k )  

        # check convergence ... 
        dx = norm(x_min - x_k) 
        if  ( dx < tol ) && 
            ( norm(c_fn(x_min)) < tol ) && 
            ( norm(h_fn(x_min)) < tol ) 
                loop = false 
        else 
            x_k = x_min 
        end 

        # update multipliers!!!  

        # first split λ and p into eq and ineq constraints 
        λ_c = λ_k[1:N_c] ;  λ_h = λ_k[N_c+1:end] 
        p_c = p_k[1:N_c] ;  p_h = p_k[N_c+1:end] 

        # deal with equality constraints first 
        if norm( c_fn(x_k) ) > tol 
            λ_c += p_c .* c_fn(x_k) 
            p_c *= γ 
        end 

        # now inequality constraints 
        h_k = h_fn( x_k ) 
        for i in eachindex( λ_h )
                
            # update λ
            λ_h[i] = max( λ_h[i] + p_h[i] * h_k[i] , 0.0 ) 

            # update p 
            if h_k[i] > 0 
                p_h[i] *= γ 
            end 

        end 

        # assign multipliers back in 
        λ_k = [ λ_c ; λ_h ] 
        p_k = [ p_c ; p_h ] 

    end 

    return x_k
end 

export min_aug_L_eq_ineq 