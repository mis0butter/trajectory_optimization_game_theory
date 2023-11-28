using trajectory_optimization_game_theory
using LinearAlgebra 
using ForwardDiff 
using GLMakie 

## ============================================ ##
# the big test ... 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# obj fn - true min at [-1, 0, 0]
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# eq constraints 
c_fn(x) = x[3] - 1          # x[3] = 1 

# ineq constraints: h_fn formulated as <= 0 
h_fn(x) = [ -x[1] + 1 ;     # x[1] >= 1  
            -x[2] - 1 ]     # x[2] >= -1  


x_k_check = min_aug_L_eq_ineq( obj_fn, c_fn, h_fn, x_0 ) 






## ============================================ ##
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

# assign  
aug_L_fn( obj_fn, c_fn, x_0, λ_0, p_0 ) 
fn(x) = aug_L_fn( obj_fn, c_fn, x, λ_0, p_0 ) 

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

println( "x min = ", x_k ) 

## ============================================ ## 

x_k_check = min_aug_L_eq( obj_fn, c_fn, x_0 ) 

## ============================================ ##
## ============================================ ##
# let's deal with inequality constraints  

# obj fn 
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# ineq constraints: h_fn formulated as <= 0 
h_fn(x) = [ -x[1] + 1 ;     # x[1] >= 1 --> -x[1] + 1 <= 0 
            -x[2] - 1 ]     # x[2] >= 1 --> -x[2] + 1 <= 0 
N_h     = length( h_fn(x_0) ) 

# lagrange multipliers and penalty parameters 
λ_0 = zeros(N_h) 
p_0 = 10.0 * ones(N_h) 
γ   = 2.0 

# define tol 
tol = 1e-6 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# assign augmented objective function 
fn(x) = aug_L_fn( obj_fn, h_fn, x, λ_0, p_0 ) 

# gradient fn and compute 
dfn  = x -> ForwardDiff.gradient( fn, x ) 
g    = dfn( x_0 ) 

## ============================================ ##

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
    dx = norm(x_min - x_k) ;    h_k = h_fn( x_k )
    if dx < tol 
        loop = false 
    else 
        x_k = x_min 
    end 

    # update constraint values 
    λ_k, p_k = update_λ_p_ineq( λ_k, p_k, h_k, γ ) 

end 

println( "x min = ", x_k ) 

## ============================================ ##

x_k_check = min_aug_L_ineq( obj_fn, h_fn, x_0 ) 

## ============================================ ##
## ============================================ ##
# let's deal with equality and inequality constraints!!! 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# obj fn 
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# eq constraints 
c_fn(x) = x[3] - 1 
N_c     = length( c_fn(x_0) ) 

# ineq constraints: h_fn formulated as <= 0 
h_fn(x) = [ -x[1] + 1 ;     # x[1] >= 1  
            -x[2] - 1 ]     # x[2] >= -1  
N_h     = length( h_fn(x_0) ) 

# put them together 
ψ_fn(x) = [ c_fn(x) ; h_fn(x) ]  
N_ψ     = length( ψ_fn(x_0) ) 

# lagrange multipliers and penalty parameters 
λ_0 = zeros(N_ψ) 
p_0 = 10.0 * ones(N_ψ) 
γ   = 2.0 

# define tol 
tol = 1e-6 

# first split λ and p into eq and ineq constraints 
λ_k = λ_0 
p_k = p_0 
λ_c = λ_k[1:N_c] 
p_c = p_k[1:N_c] 
λ_h = λ_k[N_c+1:end] 
p_h = p_k[N_c+1:end] 

# assign  
# fn(x) = aug_L_fn( obj_fn, ψ_fn, x, λ_0, p_0 )
fn(x) = aug_L_eq_ineq_fn( obj_fn, c_fn, h_fn, x, λ_c, λ_h, p_c, p_h )

# gradient fn and compute 
dfn  = x -> ForwardDiff.gradient( fn, x ) 
g    = dfn( x_0 ) 

## ============================================ ##
 
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

println( "x min = ", x_k ) 

## ============================================ ##

x_k_check = min_aug_L_eq_ineq( obj_fn, c_fn, h_fn, x_0 ) 

