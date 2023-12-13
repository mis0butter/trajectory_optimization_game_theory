# Lambert
  
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
    # rv_f_kepler = prop_kepler_tof( rv_0_lb, tof, mu ) 
    t, rv_prop  = propagate_2Body( rv_0_lb, tof, mu ) 
    rv_prop     = vv2m(rv_prop) 

    return rv_prop, Δv_vec 
end 

export prop_lambert_soln 

## ============================================ ##

"Minimizing Δv for Lambert solution using crappy grid search"
function crappy_min_lambert( 
    rv_0,           # initial position vector 
    rv_f,           # final position vector 
    mu  = 1.0,      # gravitational parameter 
    dm  = "pro",    # direction of motion 
)

    # get max period  
    T_0   = orbitPeriod( cart2kep(rv_0, mu), mu)  
    T_f   = orbitPeriod( cart2kep(rv_f, mu), mu) 
    T_max = maximum( [ T_0, T_f ] ) 

    # vary tof for lambert, get min Δv 
    Δv_norm = [] 
    T_vec   = collect( T_max / 10 : 100 : T_max ) 
    for tof_i = T_vec  
        _, Δv  = prop_lambert_soln( rv_0, rv_f, tof_i, dm, mu )
        push!( Δv_norm, norm(Δv) ) 
    end 

    i_min = get_index( Δv_norm, minimum(Δv_norm) ) 
    tof   = T_vec[i_min]

    return tof 
end 

export crappy_min_lambert 

## ============================================ ##

"Set initial guess for optimization using lambert solution"
function lambert_init_guess( 
    rv_0,           # initial position vector 
    rv_f,           # final position vector 
    N   = 20,       # number of segments 
    mu  = 1.0,      # gravitational parameter 
    dm  = "pro",    # direction of motion 
)

    tof = crappy_min_lambert( rv_0, rv_f, mu ) 

    # first compute lambert 
    _, Δv  = prop_lambert_soln( rv_0, rv_f, tof, dm, mu )

    # set initial guess 
    tof_N  = tof / N 

    Δv_vec = [ ]
    for i = 1 : N 
        push!( Δv_vec, Δv / N ) 
    end 
    Δv_vec      = vv2m( Δv_vec ) 

    return tof_N, Δv_vec
end 

export lambert_init_guess 

