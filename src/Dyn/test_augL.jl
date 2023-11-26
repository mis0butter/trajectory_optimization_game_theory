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
λ_0 = [ 0.0; 0.0 ]  
p_0 = [ 10.0 ; 10.0 ]  
γ   = 2.0 

# define tol 
tol = 1e-6 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# augmented lagrangian 
function augL_fn(x, λ, p, γ) 

    augL = obj_fn(x) + λ' * c_fn(x) + (p./2)' * c_fn(x).^2 

    return augL 
end 
augL_fn( x_0, λ_0, p_0, γ ) 

# assign  
fn(x) = augL_fn(x, λ_0, p_0, γ_0) 

# gradient fn 
dfn  = x -> ForwardDiff.gradient( fn, x ) 

# compute gradient 
g  = dfn( x_0 ) 

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
    fn(x_k) = augL_fn(x_k, λ_k, p_k, γ) 
    dfn  = x_k -> ForwardDiff.gradient( fn, x_k ) 

    # step 2: minimize unconstrained problem  
    x_min = min_bfgs( fn, dfn, x_k )  
    println( "x min = ", x_min ) 

    # step 3: check constraint function and update parameters 
    if norm( c_fn(x_min) ) > tol 
        λ_k += p_k .* c_fn(x_min) 
        p_k *= γ 
    else 
        loop = false 
    end 

    # step 4: update x 
    x_k = x_min 

end 


## ============================================ ##
# let's deal with one inequality constraint 

# want x[1] >= 0 . h_fn formulated as <= 0 
h_fn(x) = -x[1] + 1

# step 0: initialize 
λ_k = copy( λ_0 )
p_k = copy( p_0 ) 
x_k = copy( x_0 ) 
h_k = h_fn( x_k )

k = 0 ; loop = true 
while loop 

    k += 1 

    # step 1: assign augmented Lagrangian fn 
    if h_k < -λ_k / p_k 
        fn(x_k) = obj_fn(x_k) - λ_k^2 / (2*p_k) 
    else 
        fn(x_k) = obj_fn(x_k) + λ_k' * h_fn(x_k) + (p_k/2)' * h_k^2 
    end
    dfn  = x_k -> ForwardDiff.gradient( fn, x_k ) 

    # step 2: minimize unconstrained problem  
    x_min = min_bfgs( fn, dfn, x_k )  
    println( "x min = ", x_min ) 

    # step 3 check convergence ... 
    dx = norm(x_min - x_k) 
    if dx < tol 
        loop = false 
    end 

    # ... and update parameters 
    λ_k = max( λ_k + p_k' * h_k, 0.0 ) 
    h_k = h_fn( x_min ) 
    if h_k > 0 
        p_k *= γ
    end 

    # step 4: update x 
    x_k = x_min 

end 

# ----------------------- #
# plot for sanity check 

fn(x,y) = obj_fn([x,y]) 
x, y = -2:0.1:2, -2:0.1:2 
z_obj = fn.(x,y')

# plot 
fig = plot_surface( x, y, z_obj ) 
fig = plot_scatter( x_min[1], x_min[2], obj_fn(x_min), fig ) 
