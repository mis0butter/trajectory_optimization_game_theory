using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences
using LinearAlgebra 

## ============================================ ##
# init 

# ----------------------- #
# define IC and target state 

mu = 398600.4415
r  = 6378.0
kep0_P = [r+400.0, 0.0, 0*pi/180, 0.0, 0.0, 0.0]
kep0_E = [r+450.0, 0.0, 51.6*pi/180, 0.0, 0.0, 90*pi/180]
t = (0.0, 1*orbitPeriod(kep0_E, mu)) 
t_P, x_P = propagate_2Body(kep2cart(kep0_P, mu), t, mu, 1.0)
t_E, x_E = propagate_2Body(kep2cart(kep0_E, mu), t, mu, 1.0)

x_0 = x_P[1] 

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