using trajectory_optimization_game_theory 
using ForwardDiff 
using FiniteDifferences
using LinearAlgebra 

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


