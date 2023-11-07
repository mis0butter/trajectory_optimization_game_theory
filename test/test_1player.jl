using trajectory_optimization_game_theory
using LinearAlgebra, StaticArrays

# Infiltrator.clear_disabled!()
# Infiltrator.toggle_async_check(false) 

## ============================================ ##

mu = 398600.4415
r = 6378.0
kep0_p1 = [r+400.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]
kep0_p2 = [r+450.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]
t = (0.0, orbitPeriod(kep0_p2, mu))
prop1 = propagate_2Body(kep2cart(kep0_p1, mu), t, mu)
prop2 = propagate_2Body(kep2cart(kep0_p2, mu), t, mu)

x₀_P = prop1.u[1] 
x₀_E = prop2.u[1] 
xf_E = prop1.u[end] 

x0  = x₀_P 
xfₒ = xf_E 

sf  = solve_transfer(x0, 100, xfₒ, t0, mu)
fig = plot_solution!(x0, sf.xf, sf.Δτ, sf.Δv⃗, mu) 



