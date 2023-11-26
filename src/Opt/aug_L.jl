using LinearAlgebra 

## ============================================ ##

"Construct Augmented Lagrangian function"
function aug_L_fn( 
    x, 
    λ, 
    p, 
    obj_fn, 
    ψ_fn 
) 

    # augL = obj_fn(x) + λ' * ψ_fn(x) + (p./2)' * ψ_fn(x).^2 

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
    c_fn,                                           # constraint function 
    x_0,                                            # initial guess 
    λ_0     = zeros(length(c_fn(x_0))),             # initial Lagrange multiplier  
    p_0     = 10.0 * ones(length(c_fn(x_0))),       # initial penalty parameter  
    γ       = 2.0,                                  # penalty parameter 
    tol     = 1e-6, 
) 

    # step 0: initialize 
    λ_k = copy( λ_0 )
    p_k = copy( p_0 ) 
    x_k = copy( x_0 ) 

    k = 0 ; loop = true 
    while loop 

        k += 1 

        # step 1: assign augmented Lagrangian fn 
        fn(x_k) = aug_L_fn( x_k, λ_k, p_k, γ, obj_fn, c_fn ) 
        dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

        # step 2: minimize unconstrained problem  
        x_k = min_bfgs( fn, dfn, x_k )  

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



