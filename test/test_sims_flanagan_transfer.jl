using trajectory_optimization_game_theory
using LinearAlgebra, StaticArrays

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

sf  = sims_flanagan_transfer(x0, 100, x₀_P, t0, μ)
fig = plot_solution!(x0, sf.xf, sf.Δτ, sf.Δv⃗, μ) 
