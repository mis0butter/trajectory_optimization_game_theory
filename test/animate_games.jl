"""
Pursuit-Evasion Game Animation Script

This script simulates a pursuit-evasion game in orbital mechanics where:
- A pursuer (P) attempts to intercept an evader (E)
- Both vehicles follow two-body orbital dynamics
- The pursuer uses Model Predictive Control (MPC) with N-segment trajectory optimization
- An animation is generated showing the trajectories over time

The MPC approach:
1. At each time step, predict evader's future position
2. Solve for optimal pursuer trajectory to intercept
3. Execute first segment of trajectory
4. Repeat with updated positions
"""

using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences 
using LinearAlgebra 
using Optim 
using Plots 

## ====================================================================
## 1. INITIAL PARAMETERS AND SETUP
## ====================================================================

# Gravitational parameter (km^3/s^2) - Earth's standard gravitational parameter
mu = 398600.4415

# Earth radius (km)
r = 6378.0

# Initial orbital elements for Pursuer (P)
# Format: [semi-major axis (km), eccentricity, inclination (rad), 
#          RAAN (rad), argument of periapsis (rad), true anomaly (rad)]
kep0_P = [r + 400.0, 0.1, -20*pi/180, 10.0*pi/180, 20.0*pi/180, 30.0*pi/180]
rv_0_P = rv_0_P_OG = kep2cart(kep0_P, mu)  # Convert to Cartesian (position-velocity)

# Initial orbital elements for Evader (E)
kep0_E = [r + 450.0, 0.2, 10.6*pi/180, 40.0*pi/180, 0.0, 40.0*pi/180]
rv_0_E = rv_0_E_OG = kep2cart(kep0_E, mu)

# Time of flight for interception (seconds)
tof = 2000

# Propagate both vehicles forward to see their unperturbed trajectories
t_E, rv_E = propagate_2Body(rv_0_E, tof, mu, 1.0)
t_P, rv_P = propagate_2Body(rv_0_P, tof, mu, 1.0)
rv_P = vv2m(rv_P)  # Convert vector of vectors to matrix
rv_E = vv2m(rv_E)

# Plot initial orbits (optional visualization)
fig = plot_axes3d()
fig = plot_orbit(rv_P, fig)
fig = plot_orbit(rv_E, fig)

## ====================================================================
## 2. TRAJECTORY OPTIMIZATION SETUP
## ====================================================================
## Break trajectory into N segments for optimization

# Target: evader's final position
rv_f = rv_E[end, :]
rv_0 = rv_0_P

# Number of segments for trajectory optimization
N = 20
tof_N = tof / N  # Time per segment

# Solve for optimal delta-V sequence to intercept evader
# min_Δv: minimizes total delta-V magnitude
# min_Δv_dist: minimizes distance at final time (alternative method)
Δv_sol = min_Δv(rv_0, rv_f, tof, N, mu)
Δv_sol2 = min_Δv_dist(rv_0, rv_f, tof, N, mu)  # Not used, kept for reference

# Propagate optimized trajectory
t, rv_kepler = prop_kepler_tof_Nseg(rv_0, Δv_sol, N, tof / N, mu)

## ====================================================================
## 3. MODEL PREDICTIVE CONTROL (MPC) SIMULATION
## Simulate the pursuit-evasion game using MPC:
## - At each step, pursuer replans trajectory based on current evader prediction
## - Executes first portion of planned trajectory
## - Repeats with updated positions
## ====================================================================

# Reset to initial conditions
rv_0_E = copy(rv_0_E_OG)
rv_0_P = copy(rv_0_P_OG)

# History arrays to store trajectories for animation
rv_P_dtsim_hist = []  # Pursuer actual trajectory (executed segments)
t_P_dtsim_hist = []
rv_E_dtsim_hist = []  # Evader actual trajectory
t_E_dtsim_hist = []
rv_P_tof_hist = []    # Pursuer predicted full trajectory (for visualization)
t_P_tof_hist = []
rv_E_tof_hist = []    # Evader predicted full trajectory
t_E_tof_hist = []

# MPC parameters
N_prop = 4              # Number of segments to execute before replanning
dt_sim = tof_N * N_prop # Simulation time step (time between replanning)
k_sim = 0               # Simulation step counter

# MPC loop: run for 4 replanning cycles
for t_sim = dt_sim : dt_sim : dt_sim * 4
    
    k_sim += 1
    println("k_sim = ", k_sim)
    
    # ---------------------------------- 
    # Step 1: Predict evader's future trajectory
    # ---------------------------------- 

    # Propagate evader forward for full time-of-flight to predict where it will be
    t_E_tof, rv_E_tof = propagate_2Body(rv_0_E, tof, mu, 1.0)
    rv_E_tof = vv2m(rv_E_tof)
    push!(rv_E_tof_hist, rv_E_tof)  # Store for visualization
    push!(t_E_tof_hist, t_E_tof)
    
    # Extract predicted final position of evader
    rv_f_E = rv_E_tof[end, :]
    
    # ---------------------------------- 
    # Step 2: Solve for pursuer's optimal trajectory
    # ---------------------------------- 
    
    # Compute optimal delta-V sequence to intercept predicted evader position
    Δv_P = min_Δv(rv_0_P, rv_f_E, tof, N, mu)
    # Alternative: Δv_P = min_Δv_dist(rv_0_P, rv_f_E, tof, N, mu)
    
    # ---------------------------------- 
    # Step 3: Compute full predicted trajectory (for visualization)
    # ---------------------------------- 
    
    # Generate complete pursuer trajectory over full time-of-flight
    t_P_tof, rv_P_tof = prop_2Body_tof_Nseg(rv_0_P, Δv_P, N, tof_N, mu)
    push!(rv_P_tof_hist, rv_P_tof)  # Store for visualization
    push!(t_P_tof_hist, t_P_tof)
    
    # ---------------------------------- 
    # Step 4: Execute first N_prop segments of planned trajectory
    # ---------------------------------- 

    # Actually move pursuer forward by executing first few segments
    t_P_dtsim, rv_P_dtsim = prop_2Body_tof_Nseg(rv_0_P, Δv_P, N_prop, tof_N, mu)
    push!(rv_P_dtsim_hist, rv_P_dtsim)  # Store actual executed trajectory
    push!(t_P_dtsim_hist, t_P_dtsim)
    
    # Update pursuer's current state for next iteration
    rv_0_P = rv_P_dtsim[end, :]
    
    # ---------------------------------- 
    # Step 5: Update evader position
    # ---------------------------------- 

    # (Optional: apply evader maneuver here)
    # Δi = 5.0*pi/180
    # Δv_E = computeInclinationChange(rv_0_E, Δi, mu)
    # rv_0_E += [0.0, 0.0, 0.0, Δv_E[1], Δv_E[2], Δv_E[3]]
    
    # Propagate evader forward by simulation time step
    t_E_dtsim, rv_E_dtsim = propagate_2Body(rv_0_E, dt_sim, mu, 1.0)
    rv_E_dtsim = vv2m(rv_E_dtsim)
    push!(rv_E_dtsim_hist, rv_E_dtsim)  # Store actual trajectory
    push!(t_E_dtsim_hist, t_E_dtsim)
    
    # Update evader's current state for next iteration
    rv_0_E = rv_E_dtsim[end, :]
    
end

## ====================================================================
## 4. CREATE ANIMATION
## ====================================================================
## Generate animated GIF showing pursuit-evasion game evolution

# Set plot styling
plot_font = "Computer Modern"
default(
    fontfamily = plot_font,
    linewidth = 2,
    framestyle = :box,
    label = nothing,
    grid = false,
    markerstrokewidth = 0,
)

# Initialize animation
a = Animation()

# Create frame for each simulation step (starting from step 2)
for i = 2 : k_sim
    
    # ---------------------------------- 
    # Build cumulative trajectory history
    # ---------------------------------- 

    # Concatenate all executed trajectory segments up to current step
    rv_P_dtsim = rv_P_dtsim_hist[1]  # Start with first segment
    rv_E_dtsim = rv_E_dtsim_hist[1]
    for j = 2 : i
        rv_P_dtsim = [rv_P_dtsim; rv_P_dtsim_hist[j]]  # Append subsequent segments
        rv_E_dtsim = [rv_E_dtsim; rv_E_dtsim_hist[j]]
    end
    
    # Get predicted trajectories for current step (for visualization)
    rv_P_tof = rv_P_tof_hist[i]  # Pursuer's predicted full trajectory
    rv_E_tof = rv_E_tof_hist[i]  # Evader's predicted full trajectory
    
    # ---------------------------------- 
    # Create 3D plot
    # ---------------------------------- 

    lims = 1.1 .* (-r, r)  # Plot limits (slightly larger than Earth radius)
    plt = plot3d(
        1,
        # xlim = lims,  # Uncomment to set fixed axis limits
        # ylim = lims,
        # zlim = lims,
        title = "Pursuit Evasion Game",
        legend = false,
        xlabel = "x (km)",
        guidefontsize = 10,
        tickfontsize = 8,
        ylabel = "y (km)",
        zlabel = "z (km)",
        titlefont = font(16, "Computer Modern"),
        # camera = (-30, 35, 30),  # Uncomment to set camera angle
    )
    
    # Optional: Add Earth sphere
    # Plots.surface!(
    #     sphere(r, zeros(3)),
    #     alpha = 0.1
    # )
    
    # ---------------------------------- 
    # Plot Pursuer trajectories
    # ---------------------------------- 

    # Actual executed trajectory (solid blue line)
    plot3d!(
        rv_P_dtsim[:, 1], rv_P_dtsim[:, 2], rv_P_dtsim[:, 3],
        legend = false,
        color = :blue,
    )
    # Predicted full trajectory (semi-transparent cyan line)
    plot3d!(
        rv_P_tof[:, 1], rv_P_tof[:, 2], rv_P_tof[:, 3],
        linealpha = 0.5,
        color = :cyan,
    )
    # Starting position marker (x-cross)
    scatter3d!(
        [rv_P_dtsim[1, 1]], [rv_P_dtsim[1, 2]], [rv_P_dtsim[1, 3]],
        markershape = :xcross,
        color = :blue,
    )
    # Current position marker (triangle)
    scatter3d!(
        [rv_P_dtsim[end, 1]], [rv_P_dtsim[end, 2]], [rv_P_dtsim[end, 3]],
        markershape = :utriangle,
        color = :blue,
    )
    
    # ---------------------------------- 
    # Plot Evader trajectories
    # ---------------------------------- 

    # Actual executed trajectory (solid red line)
    plot3d!(
        rv_E_dtsim[:, 1], rv_E_dtsim[:, 2], rv_E_dtsim[:, 3],
        color = :red,
    )
    # Predicted full trajectory (semi-transparent orange line)
    plot3d!(
        rv_E_tof[:, 1], rv_E_tof[:, 2], rv_E_tof[:, 3],
        linealpha = 0.5,
        color = :orange,
    )
    # Starting position marker (x-cross)
    scatter3d!(
        [rv_E_dtsim[1, 1]], [rv_E_dtsim[1, 2]], [rv_E_dtsim[1, 3]],
        markershape = :xcross,
        color = :red,
    )
    # Current position marker (triangle)
    scatter3d!(
        [rv_E_dtsim[end, 1]], [rv_E_dtsim[end, 2]], [rv_E_dtsim[end, 3]],
        markershape = :utriangle,
        color = :red,
    )
    
    # Add frame to animation
    frame(a, plt)
    
    # Save individual frame as PNG
    filename_string = string("test/outputs/test_intercept_", i, ".png")
    savefig(filename_string)
    
end

# Create and display GIF animation
g = gif(a, fps = 2.0)
display(g)
