## ====================================================================

"Minimize Δv for trajectory with N segments "
function min_Δv(
    rv_0,               # initial position vector
    rv_f,               # final position vector
    tof,                # time of flight
    N      = 20,        # number of segments
    mu     = 1.0,       # gravitational parameter
    ;                   # --- keyword-only below (see note in min_Δv_dist) ---
    dm     = "pro",     # direction of motion
    Δv_max = 2.0,       # maximum per-segment Δv [km/s]
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

## ====================================================================

"""
    min_Δv_dist(rv_0, rv_f, tof, N, mu; dm, Δv_max, w_miss)

Minimize Δv for a trajectory with N segments, targeting `rv_f`.

`dm` and `Δv_max` are **keyword-only, deliberately**. They used to be positional
arguments 6 and 7, and the only production call site (`compute_players_XU`)
passed five positional arguments — so `Δv_max` silently fell back to its 2.0 km/s
default and the cap was never threaded from `params`. Making them keywords means
a caller can set `Δv_max` without also having to know about `dm` (defect D6).

`w_miss` scales the terminal-miss term against the fuel term. This weighting is
not cosmetic: the miss term is a distance in km (order 1-40 at the Lambert
initial guess) while `sum_norm_Δv` returns an energy-like Σ‖Δv‖² in (km/s)²
(order 1e-3), so at `w_miss = 1.0` the miss term outweighs fuel by roughly three
orders of magnitude and the solve is effectively pure targeting with fuel
ignored (defect D11). Default 1.0 preserves the historical behavior exactly.
"""
function min_Δv_dist(
    rv_0,               # initial position vector
    rv_f,               # final position vector
    tof,                # time of flight
    N      = 20,        # number of segments
    mu     = 1.0,       # gravitational parameter
    ;                   # --- keyword-only below ---
    dm     = "pro",     # direction of motion
    Δv_max = 2.0,       # maximum per-segment Δv [km/s]
    w_miss = 1.0,       # weight on the terminal miss term
)

    min_Δv_dist_solve(rv_0, rv_f, tof, N, mu; dm, Δv_max, w_miss).Δv_sol
end

export min_Δv_dist

"""
    min_Δv_dist_solve(...) -> (; Δv_sol, converged, g_converged, iters, t, escalated, max_Δv)

`min_Δv_dist` plus per-solve diagnostics, which is what `compute_players_XU` records so the
convergence rate over the ~450,000 solves in a sweep can be reported (Reviewer 7 #1).

**The augmented Lagrangian is bypassed unless it is needed.** The per-segment ΔV cap is
inactive by a factor of ~90 in this scenario (measured max ‖Δv‖ ≈ 0.022 km/s against a 2.0 km/s
cap, defect D6), so wrapping every solve in an outer AL loop only adds iterations. Here the
problem is solved unconstrained, the cap is then *verified*, and the AL is invoked only if it is
genuinely violated — which also keeps the code correct if someone sets a binding `Δv_max`.
`escalated` records which path was taken.

The objective is **guarded**: a large enough trial Δv makes the Kepler propagation non-finite,
and a `NaN` trips an assertion inside Optim's line search. Returning a large finite value
instead lets the line search back off.
"""
function min_Δv_dist_solve(
    rv_0, rv_f, tof,
    N       = 20,
    mu      = 1.0,
    ;
    dm      = "pro",
    Δv_max  = 2.0,
    w_miss  = 1.0,
    method  = default_method(),
    tol     = 1e-10,
    maxiter = 2000,
)

    tof_N, Δv_vec = lambert_init_guess( rv_0, rv_f, tof, N, mu, dm )
    x_0 = vec( Δv_vec )

    raw(x) = sum_norm_Δv( x, N ) +
             w_miss * miss_distance_prop_kepler_Nseg( rv_0, x, N, rv_f, tof_N, mu )
    guarded(x) = (v = raw(x); isfinite(v) ? v : 1e6)

    # nondimensionalize: decision variables are O(1e-2) while objective gradients are
    # O(1e3), since a small Δv moves the terminal position a long way over 10 segments
    s  = max( maximum(abs, x_0), 1e-6 )
    z0 = x_0 ./ s
    r  = min_optim_info( z -> guarded(s .* z), z0; method, tol, maxiter )
    x_min = s .* vec(r.x_min)

    # verify the ΔV cap, and escalate only if it actually binds
    escalated = any( >(0), constrain_Δv( x_min, N, Δv_max ) )
    if escalated
        h_fn(x) = constrain_Δv( x, N, Δv_max )
        x_min = vec( min_aug_L( guarded, reshape(x_min, N*3, 1), nothing, h_fn ) )
    end

    Δv_sol = reshape( x_min, N, 3 )

    (; Δv_sol,
       converged   = r.converged,
       g_converged = r.g_converged,
       iters       = r.iters,
       t           = r.t,
       escalated,
       max_Δv      = maximum(norm.(eachrow(Δv_sol))))
end

export min_Δv_dist_solve

## ====================================================================

"Maximize Δv for trajectory with N segments (seems to be not working)"
function max_Δv_dist(  
    rv_0,               # initial position vector 
    rv_f,               # final position vector 
    tof,                # time of flight 
    N      = 20,        # number of segments
    mu     = 1.0,       # gravitational parameter
    ;                   # --- keyword-only below (see note in min_Δv_dist) ---
    dm     = "pro",     # direction of motion
    Δv_max = 2.0,       # maximum per-segment Δv [km/s]
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

## ====================================================================

"""
Default inner solver: **LBFGS** with a backtracking line search.

Chosen over BFGS for the Gate B sweep after the λ2 = 0 decision. Measured over 60 subproblems
spanning all game steps:

    BFGS  + BackTracking   100% converged   828 iters   miss 7.5e-12   ΔV 0.0101
    LBFGS + BackTracking   100% converged    56 iters   miss 4.2e-12   ΔV 0.0160

LBFGS is ~15x cheaper and *more* accurate on terminal miss; its only cost is 35% more ΔV. With
the game payoff now separation-only (note [j]), that ΔV difference never reaches a decision — it
changes a reported metric, not play. Swap back with `min_Δv_dist_solve(...; method = BFGS(...))`
if λ2 is ever raised, because then fuel does drive the game.
"""
default_method() = LBFGS(linesearch = Optim.LineSearches.BackTracking())

"""
    min_optim_info(fn, x_0; method, tol, maxiter) -> (; x_min, converged, g_converged, iters, t)

Gradient-based unconstrained minimization, returning solver diagnostics.

**Was Nelder-Mead (defect D9).** The old body built a ForwardDiff gradient and passed it to
`optimize(fn, dfn, x_0, NelderMead())`, which ignores it — so a derivative-free simplex method
ran on a 30-dimensional problem while an analytic gradient sat unused. It could not have been
otherwise: `ForwardDiff.gradient` *threw* until D10 was fixed in `cart2kep`.

Measured on the 12 step-1 trajectory subproblems (median):

    AL + NelderMead (incumbent)   0.122 s   miss 7.0e-05 km   ΔV 0.04758
    BFGS  + BackTracking          0.154 s   miss 7.6e-12 km   ΔV 0.02926   12/12 converged
    LBFGS + BackTracking          0.016 s   miss 3.2e-12 km   ΔV 0.04169   12/12 converged
    Ipopt via JuMP @operator      0.162 s   miss 1.0e-06 km   ΔV 0.03187    0/12 converged

BFGS is the default: 38% less fuel and ~10^7 better terminal miss than the incumbent, for 26%
more wall time. `LBFGS()` is the fast alternative if sweep time ever binds.

Two details are load-bearing:

  * **`BackTracking`, not the `HagerZhang` default.** HagerZhang extrapolates, and a large enough
    trial step makes the Kepler propagation non-finite, which trips an assertion inside the line
    search (`isfinite(phi_c)`). Backtracking only ever shrinks the step. With HagerZhang only
    1/12 solves converged; with BackTracking, 12/12.
  * **The gradient must be in-place.** `Optim` wants `g!(G, x)`; handing it an out-of-place
    closure is what let the old code silently pair a gradient with a derivative-free method.

Note `g_converged` is typically false: these terminate on step/objective tolerance rather than
gradient norm, because the terminal-miss term enters as an *exact penalty* (a norm, hence
non-smooth exactly at the solution being sought). Report that honestly rather than claiming
first-order optimality.
"""
function min_optim_info(
    fn,                     # objective function
    x_0,                    # initial guess
    ;
    method  = default_method(),
    tol     = 1e-10,
    maxiter = 2000,
)

    x0v = vec(x_0)
    cfg = ForwardDiff.GradientConfig(fn, x0v)
    g!(G, x) = (ForwardDiff.gradient!(G, fn, x, cfg); G)

    t = @elapsed result = optimize(fn, g!, x0v, method,
            Optim.Options(g_tol = tol, iterations = maxiter, allow_f_increases = true))

    # preserve the caller's array shape: aug_L passes a 30x1 Matrix and then does
    # norm(x_min - x_k), which would be a DimensionMismatch against a plain Vector
    x_min = reshape(Optim.minimizer(result), size(x_0))

    (; x_min,
       converged   = Optim.converged(result),
       g_converged = Optim.g_converged(result),
       iters       = Optim.iterations(result),
       f_min       = Optim.minimum(result),
       t)
end

export min_optim_info

"Minimize function using Optim (shape-preserving; see min_optim_info for diagnostics)."
function min_optim(
    fn,                     # objective function
    x_0,                    # initial guess
    method = default_method(),
    tol    = 1e-10,
)
    min_optim_info(fn, x_0; method, tol).x_min
end

export min_optim 

## ====================================================================

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

## ====================================================================

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


## ====================================================================

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

## ====================================================================

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
