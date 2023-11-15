using trajectory_optimization_game_theory

## ============================================ ## 

r1      = [20.0e6, 20.0e6, 0]   # [m] 
r2      = [-20.0e6, 10.0e6, 0]  # [m] 
tof     = 1.0 * 86400 
mu      = 398600.4418e9         # [m^3/s^2] 
dm      = "retro" 
Dtsec   = tof 
v1, v2  = lambertbattin(r1, r2, mu, dm, tof) 










