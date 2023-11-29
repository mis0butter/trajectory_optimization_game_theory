using trajectory_optimization_game_theory
using LinearAlgebra 
using ForwardDiff 
using GLMakie 

## ============================================ ##
# example equality and inequality min test 

# initial guess 
x_0 = [ 2.0, 2.0, 2.0 ] 

# obj fn - true min at [-1, 0, 0]
obj_fn(x) = (x[1] + 1)^2 + x[2]^2  + x[3]^2 

# eq constraints 
c_fn(x) = x[3] - 1          # x[3] = 0 --> x[3] - 1 = 0 

# ineq constraints: h_fn formulated as <= 0 
h_fn(x) = [ -x[1] + 1 ;     # x[1] >= 1   --> -x[1] + 1 <= 0    * active 
            -x[2] - 1 ]     # x[2] >= -1  --> -x[2] - 1 <= 0    * inactive 

# expect min at [1, 0, 1]
x_min = min_aug_L_eq_ineq( obj_fn, c_fn, h_fn, x_0 ) 


## ============================================ ##
# lambert 

# ----------------------- #
# define IC and target state 

mu = 398600.4415
r  = 6378.0
kep0_P = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_E = [r+450.0, 0.0, 51.6*pi/180, 0.0, 0.0, 90*pi/180]
t = (0.0, 1*orbitPeriod(kep0_E, mu)) 
t_P, x_P = propagate_2Body(kep2cart(kep0_P, mu), t, mu, 1.0)
t_E, x_E = propagate_2Body(kep2cart(kep0_E, mu), t, mu, 1.0)

x₀_P = x_P[1] 
x₀_E = x_E[1] 
xf_E = x_E[end] 

# ----------------------- #
# lambert solve 

r_0 = x₀_P[1:3] 
r_f = xf_E[1:3] 

tof     = t[end]
dm      = "retro" 
Dtsec   = tof 
v_0  = lambertbattin(r_0, r_f, mu, dm, tof) 

x₀_P_lambert   = [r_0; v_0] 
t_P_lambert, x_P_lambert = propagate_2Body( x₀_P_lambert, tof, mu )

# ----------------------- #
# vector of vectors --> matrix 

x_P = mapreduce( permutedims, vcat, x_P ) 
x_E = mapreduce( permutedims, vcat, x_E ) 
x_P_lambert = mapreduce( permutedims, vcat, x_P_lambert ) 

# ----------------------- #

fig = plot_orbit( x_P ) 
fig = plot_orbit( x_E, fig ) 
fig = plot_orbit( x_P_lambert, fig ) 
fig = plot_scatter3d( r_f[1], r_f[2], r_f[3], fig ) 

## ============================================ ##
# 
# intercept problem: 
# 
# minimize 
# 		v_0 = lambertbattin(r_0, r_f, mu, dm, tof) 
# 
# subject to: 	
# 
# 		tof <= constraint 
# 
# state: 		
# 		tof, v0 

# first guess 
# x_0 = [ tof * 1.1; 
#         v_0 + 0.5*randn(3) ] 
x_0 = tof 

# obj fn 
obj_fn(x) = norm( lambertbattin( r_0, r_f, mu, dm, x[1] ) )
obj_fn(x) = norm( lambertbattin( r_0, r_f, mu, dm, x ) )
obj_fn(x_0) 

if length(x_0) == 1 
    dfn = x -> ForwardDiff.derivative( obj_fn, x ) 
else 
    dfn = x -> ForwardDiff.gradient( obj_fn, x ) 
end 
dfn(x_0) 

# # ineq constraints: h(x) <= 0  
h_fn(x) = x[1] - 86400 
h_fn(x_0) 

x_min = min_aug_L_ineq( obj_fn, h_fn, x_0 ) 
# x_min = min_bfgs( obj_fn, dfn, x_0 )

## ============================================ ##

t = collect(10 : 10 : tof * 2) 

out = [] 
for i = 1 : length(t) 
    push!( out, obj_fn(t[i]) ) 
end 
[ out_i for out_i in out ]

using Plots 
Plots.plot(t, out)


## ============================================ ##

using Optim 

# x0 = [0.0, 0.0]
f_test(x) = (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2 
x0 = [0.0, 0.0] 

result = optimize(f_test, x0, LBFGS())
x_min = result.minimizer 

result = optimize(obj_fn, x_0, LBFGS()) 
x_min = result.minimizer 


## ============================================ ##

v_0_min = x_min[2:end] 
tof_min = x_min[1] 

x₀_P_lambert   = [r_0; v_0_min] 
t_P_lambert, x_P_lambert = propagate_2Body( x₀_P_lambert, tof_min, mu )
x_P_lambert = mapreduce( permutedims, vcat, x_P_lambert ) 

fig2 = plot_orbit( x_P ) 
fig2 = plot_orbit( x_E, fig2 ) 
fig2 = plot_orbit( x_P_lambert, fig2 ) 
fig2 = plot_scatter3d( r_f[1], r_f[2], r_f[3], fig2 ) 




        




