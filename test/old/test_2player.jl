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

# Creating Initial Condition
#   Orbits mi, so basically has its orbit around Jupiter
xp = pcm2cart(propagate_PlanetaryConstantsModel(mi, t0, μ), μ) |> collect
v̂p = xp[4:6] / norm(xp[4:6])
x0 = xp + vcat(zeros(3), v̂p*v∞) 

# Solving Transfer to Europa
x0 = convert(Vector, x0) 

#   xfₒ - final state at end of time period
#   x̄f₀ - non-dimensionalized final state
xfₒ  = pcm2cart(propagate_PlanetaryConstantsModel(C.Europa, t0, μ), μ) |> collect
xfₒ  = convert(Vector, xfₒ) # QUICK FIX

# sf  = solve_transfer(x0, 100, C.Europa, t0, μ)
sf  = solve_transfer(x0, 100, xfₒ, t0, μ)
fig = plot_solution!(x0, sf.xf, sf.Δτ, sf.Δv⃗, μ) 

## ============================================ ## 


mi = C.Io 

# Creating Initial Condition
#   Orbits mi, so basically has its orbit around Jupiter
xp = pcm2cart(propagate_PlanetaryConstantsModel(mi, t0, μ), μ) |> collect
v̂p = xp[4:6]/norm(xp[4:6])
x0 = xp + vcat(zeros(3), v̂p*v∞)

# Solving Transfer to Europa
x0 = convert(Vector, x0)

# sf  = solve_transfer(x0, 100, C.Europa, t0, μ)
sf  = solve_transfer(x0, 100, xfₒ, t0, μ)
fig = plot_solution!(x0, sf.xf, sf.Δτ, sf.Δv⃗ , μ, nothing, nothing, fig)



