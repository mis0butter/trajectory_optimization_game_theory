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
r = 6378.0
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

r1 = x₀_P[1:3] 
r2 = xf_E[1:3] 

tof     = t[end]
dm      = "retro" 
Dtsec   = tof 
v1, v2  = lambertbattin(r1, r2, mu, dm, tof) 

x₀_P_lambert   = [r1; v1] 
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

## ============================================ ##
# 
# intercept problem: 
# 
# minimize 
# 		delta_v_1 
# 
# subject to: 	
# 		tof = min( tof_pro, tof_retro ) 
# 		v0  = min( v0_pro, v0_retro ) 
# 
# 		tof_pro,   v0_pro   = lambert( r_init, r_target, 'pro' )
# 		tof_retro, v0_retro = lambert( r_init, r_target, 'pro' )
# 
# 		delta_v <= constraint 
# 
# state: 		
# 		tof, v0 

tof = 1.0 
v0  = [ 1.0, 0.0, 0.0 ] 





