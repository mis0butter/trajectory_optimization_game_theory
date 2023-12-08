# =====================================================================
# === Lambert
  
include("lambertbattin.jl")
include("seebatt.jl")
include("seebattk.jl") 

# Exports
export lambertbattin, seebatt, seebattk 

## ============================================ ##

"Propagate Keplerian orbit using Lambert solution"
function prop_lambert_soln(  
    rv_0,           # initial position vector 
    rv_f,           # final position vector 
    tof,            # time of flight 
    dm  = "pro",    # direction of motion 
    mu  = 1.0,      # gravitational parameter 
)

    r_0 = rv_0[1:3] ;   v_0 = rv_0[4:6] 
    r_f = rv_f[1:3] ;   v_f = rv_f[4:6] 

    # lambert soln 
    v_0_lb, v_f_lb = lambertbattin(r_0, r_f, mu, dm, tof) 
    rv_0_lb        = [r_0 ; v_0_lb] 

    # get delta v 
    Δv_vec = v_0_lb - v_0 

    # propagate 
    rv_f_kepler = prop_kepler_tof( rv_0_lb, tof, mu ) 
    t, rv_prop  = propagate_2Body( rv_0_lb, tof, mu ) 
    rv_prop     = vv2m(rv_prop) 

    return rv_prop 
end 

export prop_lambert_soln 



