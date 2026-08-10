using DifferentialEquations
using LinearAlgebra
using Debugger


## ====================================================================

using Infiltrator

function prop_chosen_rv(rv_E, rv_P, params)

    # @infiltrate 

    t_E, rv_E_hist = prop_kepler_tof_Nseg(rv_E, zeros(params.n_seg_horizon, 3), params.n_seg_horizon, params.t_horizon / params.n_seg_horizon, params.mu)
    t_P, rv_P_hist = prop_kepler_tof_Nseg(rv_P, zeros(params.n_seg_horizon, 3), params.n_seg_horizon, params.t_horizon / params.n_seg_horizon, params.mu)

    # t_E, rv_E_hist = propagate_2Body(rv_E, params.t_horizon, params.mu, 1.0) 
    # t_P, rv_P_hist = propagate_2Body(rv_P, params.t_horizon, params.mu, 1.0) 
    # rv_P_hist = vv2m(rv_P_hist) 
    # rv_E_hist = vv2m(rv_E_hist) 

    return t_E, rv_E_hist, t_P, rv_P_hist
end

export prop_chosen_rv


## ====================================================================

"propagate orbit based on given initial conditions, time, and gravitational parameter"
function propagate_2Body(x0, t, mu=1.0, dt=nothing)

    prob = ODEProblem(eom_2Body!, x0, t, mu)

    if isnothing(dt)
        sol = solve(prob)
    else
        sol = solve(prob, saveat=dt)
    end

    t = sol.t
    x = sol.u

    return t, x
end

export propagate_2Body

## ====================================================================

export eom_2Body!
function eom_2Body!(dx, x, mu, t)
    x1 = x[1]
    x2 = x[2]
    x3 = x[3]

    mu_div_r3 = -mu / sqrt(x1^2 + x2^2 + x3^2)^3

    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] = mu_div_r3 * x1
    dx[5] = mu_div_r3 * x2
    dx[6] = mu_div_r3 * x3
end

## ====================================================================

export kep2cart
function kep2cart(kep, mu)

    sma = kep[1]
    ecc = kep[2]
    inc = kep[3]
    raan = kep[4]
    omega = kep[5]
    theta = kep[6]

    rcoeff = (sma * (1.0 - ecc^2)) / (1.0 + ecc * cos(theta))
    rx_pqw = rcoeff * cos(theta)
    ry_pqw = rcoeff * sin(theta)
    r_pqw = [rx_pqw, ry_pqw, 0]

    vcoeff = sqrt(mu / (sma * (1.0 - ecc^2)))
    vx_pqw = vcoeff * (-sin(theta))
    vy_pqw = vcoeff * (ecc + cos(theta))
    v_pqw = [vx_pqw, vy_pqw, 0]

    rot_PQW_to_IJK = R3(-raan) * R1(-inc) * R3(-omega)
    r = rot_PQW_to_IJK * r_pqw
    v = rot_PQW_to_IJK * v_pqw

    cart = [r; v]
    return cart
end

## ====================================================================

function R3(angle)
    R = [cos(angle) sin(angle) 0.0;
        -sin(angle) cos(angle) 0.0;
        0.0 0.0 1.0]
    return R
end

## ====================================================================

function R1(angle)
    R = [1.0 0.0 0.0;
        0.0 cos(angle) sin(angle);
        0.0 -sin(angle) cos(angle)]
    return R
end

## ====================================================================

export cart2kep
function cart2kep(rv, mu)

    r = rv[1:3]
    v = rv[4:6]

    pi2 = 2.0 * π

    # position and velocity magnitude
    rmag = norm(r)
    vmag = norm(v)

    # position unit vector
    rhat = r / rmag

    # angular momentum vectors
    hv = cross(r, v)
    hhat = hv / norm(hv)

    # eccentricity vector
    vtmp = v / mu
    ecc = cross(vtmp, hv) - rhat

    # semimajor axis
    sma = 1.0 / (2.0 / rmag - vmag^2 / mu)
    p = hhat[1] / (1.0 + hhat[3])
    q = -hhat[2] / (1.0 + hhat[3])

    const1 = 1.0 / (1.0 + p^2 + q^2)

    fhat = [const1 * (1.0 - p^2 + q^2),
        const1 * 2.0 * p * q,
        -const1 * 2.0 * p]

    ghat = [const1 * 2.0 * p * q,
        const1 * (1.0 + p^2 - q^2),
        const1 * 2.0 * q]

    h = dot(ecc, ghat)
    xk = dot(ecc, fhat)
    x1 = dot(r, fhat)
    y1 = dot(r, ghat)

    # orbital eccentricity
    eccm = sqrt(h^2 + xk^2)

    # orbital inclination
    inc = 2.0 * atan(sqrt(p^2 + q^2))

    # true longitude
    xlambdat = atan(y1, x1)

    # check for equatorial orbit
    raan = inc > 1e-8 ? atan(p, q) : 0.0

    # check for circular orbit
    argper = eccm > 1e-8 ? mod(atan(h, xk) - raan, pi2) : 0.0

    # true anomaly
    tanom = mod(xlambdat - raan - argper, pi2)

    # load orbital element vector
    oe = zeros(6)
    oe[1] = sma
    oe[2] = eccm
    oe[3] = inc
    oe[4] = raan
    oe[5] = argper
    oe[6] = tanom

    return oe
end


## ====================================================================

# export cart2kep
# function cart2kep(cart, mu, tol=1e-20)

#     # get position and velocity components
#     r = cart[1:3]
#     v = cart[4:6]

#     # orbit energy calc
#     h = cross(r, v) # specific angular momentum
#     # magR = norm(r)
#     magR = sqrt(r[1]^2 + r[2]^2 + r[3]^2) 
#     # magV = norm(v)
#     magV = sqrt(v[1]^2 + v[2]^2 + v[3]^2) 
#     energy = 0.5 * (magV ^ 2) - (mu / magR) # vis viva equation

#     # eccentricity vector
#     ecc = (cross(v, h) - mu * (r / magR)) / mu
#     magEcc = norm(ecc) 

#     # semi-latus rectum 
#     h_norm2 = h[1]^2 + h[2]^2 + h[3]^2 
#     p = h_norm2 / mu 

#     # semi-major axis 
#     sma = p / (1 - magEcc^2) 

#     theta = acos(dot(r, ecc) / (magR * magEcc))
#     if (dot(r, v) < 0.0)
#         theta = 2.0 * pi - theta
#     end
#     inc = acos(h[3] / norm(h))

#     K = [0.0; 0.0; 1.0]
#     n = cross(K, h)
#     normN = norm(n)

#     # get right ascension
#     raan = acos(n[1] / normN)
#     if n[2] < 0.0
#         raan = 2.0 * pi - raan
#     end

#     # get argument of periapsis 
#     out = dot(n, ecc) / (normN * magEcc) 
#     if abs(out) > 1.0
#         println( "out > 1.0" )
#         out = sign(out) 
#     end 
#     omega = acos( out )
#     if ecc[3] < 0.0
#         omega = 2.0 * pi - omega
#     end

#     # singularity checks
#     I = [1.0; 0.0; 0.0]

#     if magEcc < tol && inc < tol
#         raan = 0.0
#         omega = 0.0

#         # set true longitude of periapsis as theta
#         theta = acos(dot(r, I) / magR)
#         if r[2] < 0.0
#             theta = 2.0 * pi - theta
#         end
#     elseif magEcc < tol
#         omega = 0.0

#         # set argument of latitude as theta
#         theta = acos(dot(n, r) / (normN * magR))
#         if r[3] < 0.0
#             theta = 2.0 * pi - theta
#         end
#     elseif inc < tol
#         raan = 0.0

#         # set longitude of periapsis as omega
#         omega = acos(dot(ecc, I) / magEcc)
#         if ecc[2] < 0.0
#             omega = 2.0 * pi - omega
#         end
#     end

#     # output theta here refers to true anomaly!
#     kepState = [sma; magEcc; inc; raan; omega; theta]

#     return kepState
# end

## ====================================================================

export orbitPeriod
function orbitPeriod(kep, mu)
    T = 2.0 * pi * sqrt(kep[1]^3 / mu)
    return T
end

## ====================================================================

"Non-dimensionalize position and velocity vectors"
function nondim_rv(
    r_vec,      # position vector  
    v_vec,      # velocity vector 
    mu,         # gravitational parameter 
    R,          # Earth radius 
)

    # Distance unit DU is defined by the Earth radius 
    DU = R

    # Time unit TU is defined by Earth mu
    TU = sqrt(DU^3 / mu)

    # Converting Units
    r̄_vec = r_vec / DU
    v̄_vec = v_vec / (DU / TU)

    # Outputting
    return r̄_vec, v̄_vec, DU, TU
end

export nondim_rv

