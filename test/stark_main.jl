using trajectory_optimization_game_theory 
using ForwardDiff 
using LinearAlgebra 


## ============================================ ##
## routine for FGStark_oneStep_propagator in Julia 

# declare variables 

# X0: Initial Conditions, acc: Stark acceleration (3D), 
X0   = [ 0.2e0  0.e0  0.e0  0.e0  3.0e0  0.e0  0.e0 ] ; 
acc  = 0.1 * [ 1.e-1  1.e-1  1.e-1 ] ; 
# dtau = 1e-5 ; 

# get OEs 
mu  = 1 ; 
rv0 = X0[1:6] ; 
oe0 = cart2kep( rv0, mu ) ;  
a = oe0[1] ; e = oe0[2] ; 

# Sundman transformation parameters (dt = cSun * r**alphaSun * dtau)
cSun     = 1 ; 
alphaSun = 2 ; 

# Select a case for the Sundman transformation: 
#   0: alpha = 0, tau = time
#   1: alpha = 1, tau = eccentric anomaly
#   2: alpha = 2, tau = true anomaly
#   3: alpha = 3/2, tau = intermediate anomaly
#   4: alpha = alphaSun, generic Sundman transformation: user provided alpha
caseSun = 0 ; 

# tau period 
tau_p, alphaSun = taup_caseSun( caseSun, a, e, mu, cSun ) ; 
tau_p = 0.05 * tau_p ; 

# dtau: delta tau for integration for N steps 
# dtau = 6.283185307179589e-1 ; 
N    = 50 ; 
dtau = tau_p / N ; 

# order: order of the TS integration 
order = 8 ; 

# Xf: Integrated vector 
Xf = zeros(1,7) ; 

# Accuracy measurements 
accuracyOut = zeros(1,4) ; 

# Indicates whether the user wants a Stark or Kepler propagation 
caseStark = true ; 

# Indicates whether the user wants to perform divergence checks (expensive)
divSafeguard = true ; 

# Indicates possible divergence of the series. 
# If statusFlag = -1, the series likely diverged
statusFlag = 0 ; 

# Indicates which accuracy measurements are to be computed 
# (1 = Hamiltonian, 2 = delR, 3 = delV, 4 = delT)
accuracyFlag = [ true  true  true  true ] ; 


## ============================================ ##
## begin subroutine 

println("*** Starting the F&G Stark series propagation ***")
println("")

println("X0: ") ; println(X0)
println("")

println("Case: ", caseSun, " cSun: ", cSun, " alphaSun: ", alphaSun, " order: ", order)
println("")

# propagate N dtau steps 
X_hist, accuracy_hist, status_hist, acc_hist = FGStark_Nstep_propagator( X0, acc, dtau, N, order,  
        cSun, caseSun, alphaSun, caseStark, divSafeguard, 
        accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
# [ X_hist, accuracy_hist, status_hist, acc_hist ] = FGStark_Ntaup_propagator( ... 
#         X0, acc, dtau, N, order, mu, ... 
#         cSun, caseSun, alphaSun, caseStark, divSafeguard, ... 
#         accuracyFlag, Xf, accuracyOut, statusFlag ) ; 
    





## ============================================ ##
## ============================================ ##





