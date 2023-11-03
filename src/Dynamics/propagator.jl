using Revise
using DifferentialEquations
using LinearAlgebra 

## ============================================ ##

# propagate orbit based on given initial conditions, time, and gravitational parameter
export propagate_2Body
function propagate_2Body(x0, t, mu)
    prob = ODEProblem(eom_2Body!, x0, t, mu)
    return solve(prob)
end

## ============================================ ##

export eom_2Body! 
function eom_2Body!(dx, x, mu, t)
    x1 = x[1]
    x2 = x[2] 
    x3 = x[3]

    mu_div_r3 = -mu/sqrt(x1^2 + x2^2 + x3^2) ^ 3

    dx[1] = x[4]
    dx[2] = x[5]
    dx[3] = x[6]
    dx[4] = mu_div_r3 * x1
    dx[5] = mu_div_r3 * x2
    dx[6] = mu_div_r3 * x3
end

## ============================================ ##

export kep2cart 
function kep2cart(kep, mu)

    sma    = kep[1]
    ecc    = kep[2]
    inc    = kep[3]
    raan   = kep[4]
    omega  = kep[5]
    theta  = kep[6]

    rcoeff = (sma * (1.0 - ecc^2)) / (1.0 + ecc * cos(theta))
    rx_pqw = rcoeff * cos(theta)
    ry_pqw = rcoeff * sin(theta)
    r_pqw  = [rx_pqw, ry_pqw, 0]

    vcoeff = sqrt(mu / (sma * (1.0 - ecc^2))) 
    vx_pqw = vcoeff * (-sin(theta))
    vy_pqw = vcoeff * (ecc + cos(theta))
    v_pqw  = [vx_pqw, vy_pqw, 0]

    rot_PQW_to_IJK = R3(-raan) * R1(-inc) * R3(-omega)
    r = rot_PQW_to_IJK * r_pqw
    v = rot_PQW_to_IJK * v_pqw

    cart = [r; v]
    return cart
end

## ============================================ ##

function R3(angle)
    R = [cos(angle) sin(angle) 0.0;
         -sin(angle) cos(angle) 0.0;
         0.0 0.0 1.0]
    return R
end

## ============================================ ##

function R1(angle)
    R = [1.0 0.0 0.0;
         0.0 cos(angle) sin(angle);
         0.0 -sin(angle) cos(angle)]
    return R
end

# function cart2kep(cart, mu)

#     return kep
# end

## ============================================ ##

export orbitPeriod 
function orbitPeriod(kep, mu)
    T = 2.0 * pi * sqrt(kep[1]^3 / mu)
    return T
end

