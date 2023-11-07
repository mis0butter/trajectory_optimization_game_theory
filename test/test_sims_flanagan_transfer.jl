using trajectory_optimization_game_theory
using LinearAlgebra, StaticArrays


## ============================================ ##


mu = 398600.4415
r = 6378.0
kep0_p1 = [r+400.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]
# kep0_p2 = [42164, 0.0, 0.0, 0.0, 0.0, 0.0]
# kep0_p2 = [r+450.0, 0.0, 31.6*pi/180, 0.0, 0.0, pi]
# kep0_p2 = [r+400.0, 0.0, 51.6*pi/180, 0.0, 0.0, pi]
kep0_p2 = [r+450.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]

t = (0.0, orbitPeriod(kep0_p2, mu))
prop1 = propagate_2Body(kep2cart(kep0_p1, mu), t, mu)
prop2 = propagate_2Body(kep2cart(kep0_p2, mu), t, mu) 

x₀_P = prop1.u[1] 
x₀_E = prop2.u[1] 
xf_E = prop1.u[end] 

sf  = sims_flanagan_transfer(x₀_P, x₀_E, 100, t0, mu) 
fig = plot_sims_flanagan!(x0, sf.xf, sf.Δτ, sf.Δv⃗, mu) 


## ============================================ ##

# Infiltrator.clear_disabled!()
# Infiltrator.toggle_async_check(false) 

## ============================================ ##

# Importing Project Constants
C  = import_constants()

# Setting Initial Variables
μ  = C.μ
t0 = 1.5*86400.0
v∞ = 3.65
mi = C.Ganymede 
xp = pcm2cart(propagate_PlanetaryConstantsModel(mi, t0, μ), μ) |> collect
v̂p = xp[4:6] / norm(xp[4:6])
x₀_E = xp + vcat(zeros(3), v̂p*v∞) 

# Solving Transfer 
x₀_E = convert(Vector, x₀_E)
mi = C.Europa 
xp = pcm2cart(propagate_PlanetaryConstantsModel(mi, t0, μ), μ) |> collect
v̂p = xp[4:6] / norm(xp[4:6])

x₀_P = xp + vcat(zeros(3), v̂p*v∞) 
x₀_P = convert(Vector, x₀_P)

