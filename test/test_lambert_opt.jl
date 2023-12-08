using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 


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

## ============================================ ##
# plot 

fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 

# ----------------------- #
# lambert solution 

dm  = "pro" 

tof = t[end] / 2 
r1  = x0_P[1:3] 
r2  = xf_E[1:3] 
v1, v2 = lambertbattin(r1, r2, mu, dm, tof) 

x0_lambert  = [r1 ; v1] 

Δv_vec = v1 - x0_P[4:6] 

# test function 
rv_f_kepler = prop_kepler_tof( x0_lambert, tof, mu ) 
t, rv_prop  = propagate_2Body( x0_lambert, tof, mu ) 
rv_prop     = vv2m(rv_prop) 

println( "err norm = ", norm( rv_f_kepler[1:3] - r2 ) ) 

fig = plot_orbit( rv_prop, fig ) 
fig = plot_scatter3d( xf_E[1], xf_E[2], xf_E[3], fig ) 
fig = plot_scatter3d( rv_f_kepler[1], rv_f_kepler[2], rv_f_kepler[3], fig, :utriangle, :green ) 


## ============================================ ##
# can I use Optim? ... looks like YES 

x_0 = Δv_vec * 1.1 

# optimization 
fn(x)   = miss_distance_prop_kepler( rv_0, x, rv_f, tof, mu ) 
od      = OnceDifferentiable( fn, x_0 ; autodiff = :forward ) 
result  = optimize( od, x_0, NelderMead() ) 
x_min   = result.minimizer 

fn(x_min) 

x_min = min_optim( fn, x_0 ) 

## ============================================ ##
# ok now let's try minimizing a different objective function 

# define state vector 
x_0 = [ tof ; Δv_vec ] 

# want to minimize delta v - define objective function 
obj_fn(x) = sum( abs.( x[2:4] ) )  
obj_fn(x_0) 

# want equality constraint - miss distance = 0  
c_fn(x) = miss_distance_prop_kepler( rv_0, x[2:4], rv_f, x[1], mu ) 
c_fn(x_0) 

# x_min = min_optim( obj_fn, x_0 ) 

# equality-constrained 
x_min = min_aug_L( obj_fn, x_0, c_fn )

# extract 
tof_min = x_min[1] 
Δv_min  = x_min[2:4] 

## ============================================ ##
# propagate optimized solution to test 

# plot 
fig = plot_axes3d()
fig = plot_orbit( x_P, fig ) 
fig = plot_orbit( x_E, fig ) 

# propagate with optimized solution 
x0_min   = rv_0 + [ zeros(3) ; Δv_min ] 
rv_f_min = prop_kepler_tof( x0_min, tof_min, mu ) 

# print error stats 
println( "opt err norm = ", norm( rv_f_min[1:3] - rv_f ) ) 
println( "lambert err norm = ", norm( rv_f_kepler[1:3] - rv_f ) ) 

# test function 
t, rv_prop  = propagate_2Body( x0_min, tof_min, mu ) 
rv_prop     = vv2m(rv_prop) 

fig = plot_orbit( rv_prop, fig ) 
fig = plot_scatter3d( xf_E[1], xf_E[2], xf_E[3], fig ) 
fig = plot_scatter3d( rv_f_min[1], rv_f_min[2], rv_f_min[3], fig, :utriangle, :green ) 














## ============================================ ##


# want inequality tof constraint - less than 1 day 
h_fn(x) = x[1] - 86400 
h_fn(x_0) 

# want equality constraint - miss distance = 0  
c_fn(x) = miss_distance_prop_kepler( rv_0, x[2:4], rv_f, x[1], mu ) 
c_fn(x_0) 

# equality and inequality-constrained 
x_min = min_aug_L( obj_fn, x_0, c_fn, h_fn )



















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

