using trajectory_optimization_game_theory
using LinearAlgebra 
using ForwardDiff 
using GLMakie 

## ============================================ ##
# ok ... let's try equality constrained optimization 

# obj fn 
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# constraint fn 
c_fn(x) = [ x[2] - 1 ; 
            x[3] - 1 ] 

# lagrange multipliers and penalty parameters 
λ_0 = [ 0.0, 0.0 ]   
p_0 = [ 10.0, 10.0 ]   
γ   = 2.0 

# define tol 
tol = 1e-6 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# augmented lagrangian 
function augL_fn( x, λ, p, γ, obj_fn, ψ_fn ) 

    # augL = obj_fn(x) + λ' * ψ_fn(x) + (p./2)' * ψ_fn(x).^2 

    augL = obj_fn(x) 
    for i in eachindex(λ)
        augL += λ[i] * ψ_fn(x)[i] + (p[i]/2) * ψ_fn(x)[i]^2 
    end 

    return augL 
end 
augL_fn( x_0, λ_0, p_0, γ, obj_fn, c_fn ) 

# assign  
fn(x) = augL_fn( x, λ_0, p_0, γ, obj_fn, c_fn ) 

# gradient fn and compute 
dfn  = x -> ForwardDiff.gradient( fn, x ) 
g    = dfn( x_0 ) 

## ============================================ ##
# augmented Lagrangian method (equality-constrained) 

# step 0: initialize 
λ_k = copy( λ_0 )
p_k = copy( p_0 ) 
x_k = copy( x_0 ) 

k = 0 ; loop = true 
while loop 

    k += 1 

    # step 1: assign augmented Lagrangian fn 
    fn(x_k) = augL_fn( x_k, λ_k, p_k, γ, obj_fn, c_fn ) 
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

println( "x min = ", x_k ) 

## ============================================ ##
## ============================================ ##
# let's deal with inequality constraints  

# obj fn 
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# want: 
#   x[1] >= 1 --> -x[1] + 1 <= 0 
#   x[2] >= 1 --> -x[2] + 1 <= 0 
# h_fn formulated as <= 0 
h_fn(x) = [ -x[1] + 1 ; 
            -x[2] + 1 ] 
# h_fn(x) = -x[1] + 1  
ψ_fn    = h_fn 
N_h     = length( h_fn(x_0) ) 

# lagrange multipliers and penalty parameters 
λ_0 = zeros(N_h) 
p_0 = 10.0 * ones(N_h) 
γ   = 2.0 

# define tol 
tol = 1e-6 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# assign  
fn(x) = augL_fn( x, λ_0, p_0, γ, obj_fn, ψ_fn ) 

# gradient fn and compute 
dfn  = x -> ForwardDiff.gradient( fn, x ) 
g    = dfn( x_0 ) 

## ============================================ ##

# step 0: initialize 
λ_k = copy( λ_0 )
p_k = copy( p_0 ) 
x_k = copy( x_0 ) 
ψ_k = ψ_fn( x_k )

k = 0 ; loop = true 
while loop 

    k += 1 
    println( "k = ", k ) 

    # step 1: assign augmented Lagrangian fn 
    fn(x_k) = augL_fn( x_k, λ_k, p_k, γ, obj_fn, ψ_fn ) 
    dfn     = x_k -> ForwardDiff.gradient( fn, x_k ) 

    # step 2: minimize unconstrained problem  
    x_min = min_bfgs( fn, dfn, x_k )  

    # step 3 check convergence ... 
    dx = norm(x_min - x_k) 
    if dx < tol 
        loop = false 
    end 

    # update constraint values 
    ψ_k = ψ_fn( x_min )

    # ... and update parameters 
    for i in eachindex(λ_k)

        # update λ
        λ_k[i] = max( λ_k[i] + p_k[i] * ψ_k[i] , 0.0 ) 

        # update p 
        if ψ_k[i] > 0 
            p_k[i] *= γ 
        end 

    end 

    # step 4: update x 
    x_k = x_min 

end 

## ============================================ ##




