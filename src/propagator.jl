using Revise
using DifferentialEquations
using Plots

# propagate orbit based on given initial conditions, time, and gravitational parameter
export propagate_2Body
function propagate_2Body(x0, t, mu)
    prob = ODEProblem(eom_2Body!, x0, t, mu)
    return solve(prob)
end

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

function R3(angle)
    R = [cos(angle) sin(angle) 0.0;
         -sin(angle) cos(angle) 0.0;
         0.0 0.0 1.0]
    return R
end

function R1(angle)
    R = [1.0 0.0 0.0;
         0.0 cos(angle) sin(angle);
         0.0 -sin(angle) cos(angle)]
    return R
end

# function cart2kep(cart, mu)

#     return kep
# end

function orbitPeriod(kep, mu)
    T = 2.0 * pi * sqrt(kep[1]^3 / mu)
    return T
end

mu = 398600.4415
r = 6378.0
kep0_p1 = [r+400.0, 0.0, 51.6*pi/180, 0.0, 0.0, 0.0]
kep0_p2 = [42164, 0.0, 0.0, 0.0, 0.0, 0.0]
t = (0.0, orbitPeriod(kep0_p2, mu))
prop1 = propagate_2Body(kep2cart(kep0_p1, mu), t, mu)
prop2 = propagate_2Body(kep2cart(kep0_p2, mu), t, mu)

xVals_p1 = []
yVals_p1 = []
zVals_p1 = []
xVals_p2 = []
yVals_p2 = []
zVals_p2 = []
for i in 1:length(prob.u)
    push!(xVals_p1, prop1.u[i][1])
    push!(yVals_p1, prop1.u[i][2])
    push!(zVals_p1, prop1.u[i][3])
    push!(xVals_p2, prop2.u[i][1])
    push!(yVals_p2, prop2.u[i][2])
    push!(zVals_p2, prop2.u[i][3])
end

scatter(xVals_p1, yVals_p1, zVals_p1, label="p1")
scatter!(xVals_p2, yVals_p2, zVals_p2, label="p2")