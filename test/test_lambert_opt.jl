using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 


## ============================================ ##
# init 

mu = 398600.4415
r  = 6378.0
kep0_P = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_E = [r+450.0, 0.0, 10.6*pi/180, 0.0, 0.0, 50.0*pi/180]
t = (0.0, 1*orbitPeriod(kep0_E, mu)) 
t_P, x_P = propagate_2Body(kep2cart(kep0_P, mu), t, mu, 1.0)
t_E, x_E = propagate_2Body(kep2cart(kep0_E, mu), t, mu, 1.0)

x0_P = x0_P_OG = x_P[1] 
x0_E = x0_E_OG = x_E[1] 
xf_E = xf_E_OG = x_E[end] 

x_P = vv2m(x_P) 
x_E = vv2m(x_E) 

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 


## ============================================ ##
# lambert solution 

dm = "pro" 

tof = t[end] / 2 

r1 = x0_P[1:3] 
r2 = xf_E[1:3] 
v1, v2 = lambertbattin(r1, r2, mu, dm, tof) 

x0_lambert  = [r1 ; v1] 

Δv_vec = v1 - x0_P[4:6] 

# test function 
rv_f_kepler = prop_kepler_tof( x0_lambert, tof, mu ) 
t, rv_prop  = propagate_2Body( x0_lambert, tof, mu ) 
rv_prop     = vv2m(rv_prop) 

fig = plot_orbit( rv_prop, fig ) 
fig = plot_scatter3d( xf_E[1], xf_E[2], xf_E[3], fig ) 
fig = plot_scatter3d( rv_f_kepler[1], rv_f_kepler[2], rv_f_kepler[3], fig, :utriangle, :green ) 


## ============================================ ##
# test optimizing for initial delta v 

rv_0  = copy(x0_P) 
rv_f  = copy( r2 ) 
Δrv_f = miss_distance_prop_kepler( rv_0, Δv_vec, rv_f, tof, mu ) 

fn(x) = miss_distance_prop_kepler( rv_0, x, rv_f, tof, mu ) 
fn(Δv_vec*0.99) 

fn_Δv(x) = fn(Δv_vec * x)

x   = collect( 0.99 : 0.0001 : 1.01 ) 
out = fn_Δv.(x) 

idx_min = get_index( out, minimum(out) ) 
x[idx_min]  


## ============================================ ##

dfn     = x -> ForwardDiff.gradient( fn, x ) 
dfn_fdm = x -> grad( central_fdm(5, 1), fn, x )[1] 

dfn(Δv_vec)
dfn_fdm(Δv_vec) 

x_min_fd   = min_bfgs( fn, dfn_fdm, Δv_vec )  


## ============================================ ##
# can I use Optim? ... looks like no 

using Optim 

x0 = Δv_vec * 1.1 

# optimization 
fn(x)   = miss_distance_prop_kepler( rv_0, x, rv_f, tof, mu ) 
od      = OnceDifferentiable( fn, x0 ; autodiff = :forward ) 
result  = optimize( od, x0, NelderMead() ) 
x_min   = result.minimizer 

fn(x_min)


















## ============================================ ##
# init 

r_0, r_f, v_0, v_f, rv_lambert, Δv_vec, tof, N, mu = lambert_IC() 
# rv_0  = [ r_0 ; v_0 ]
rv_0  = [ r_0 ; 0*v_0 ] 
rv_f  = [ r_f ; v_f ] 

## ============================================ ##

# can I take gradient of lambert? 

# fn(tof) = lambertbattin( r1, r2, mu, dm, tof )
# dfn     = tof --> ForwardDiff.derivative( fn, tof ) 

fn(x) = seebatt(x)
fn(rand(1)[1]) 
dfn = x -> ForwardDiff.derivative( fn, x )  
dfn( rand(1)[1] )

fn(x) = seebattk(x)
fn(rand(1)[1]) 
dfn = x -> ForwardDiff.derivative( fn, x )  
dfn( rand(1)[1] )

fn(x) = lambertbattin( r1, r2, mu, dm, x )
fn(rand(1)[1]) 
dfn = x -> ForwardDiff.derivative( fn, x )  
dfn( rand(1)[1] )

## ============================================ ##


