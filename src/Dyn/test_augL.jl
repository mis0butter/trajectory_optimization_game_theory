using trajectory_optimization_game_theory
using LinearAlgebra 
using ForwardDiff 

## ============================================ ##
# ok ... let's try equality constrained optimization 

# obj fn 
obj_fn(x) = x[1]^2 + x[2]^2  

# constraint fn 
c_fn(x) = x[2] - 1 

# lagrange multipliers and penalty parameters 
λ_0 = 0.0 
p_0 = 10.0 
γ_0 = 2.0 

# define tol 
tol = 1e-6 

# initial guess 
x_0 = [ 2.0, 2.0 ] 

# augmented lagrangian 
function augL_fn(x, λ, p, γ) 

    return obj_fn(x) + λ' * c_fn(x) + (p./2)' * c_fn(x)^2 

end 
augL_fn( x_0, λ_0, p_0, γ_0 ) 

# assign  
fn(x) = augL_fn(x, λ, p, γ) 

# gradient fn 
dfn  = x -> ForwardDiff.gradient( fn, x ) 

# compute gradient 
g  = dfn( x_0 ) 

## ============================================ ##
# augmented Lagrangian method (equality-constrained) 

# step 0: initialize 
λ_k = copy( λ_0 )
p_k = copy( p_0 ) 
γ_k = copy( γ_0 )  
x_k = copy( x_0 ) 

k = 0 ; loop = true 
while loop 

    k += 1 

    # step 1: assign augmented Lagrangian fn 
    fn(x_k) = augL_fn(x_k, λ_k, p_k, γ_k) 
    dfn  = x_k -> ForwardDiff.gradient( fn, x_k ) 

    # step 2: minimize unconstrained problem  
    x_min = min_bfgs( fn, dfn, x_k )  
    println( "x min = ", x_min ) 

    # step 3: check constraint function and update parameters 
    if norm( c_fn(x_min) ) > tol 
        # println( "constraint not satisfied" ) 
        λ_k += p_k' * c_fn(x_min) 
        p_k *= γ_k 
    else 
        # println( "constraint satisfied" ) 
        loop = false 
    end 

    # step 4: update x 
    x_k = x_min 

end 


