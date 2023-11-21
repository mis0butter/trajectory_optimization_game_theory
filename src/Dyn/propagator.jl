using DifferentialEquations

## ============================================ ##

"propagate orbit based on given initial conditions, time, and gravitational parameter" 
function propagate_2Body(x0, t, mu = 1.0, dt = nothing)

    prob = ODEProblem(eom_2Body!, x0, t, mu)

    if isnothing(dt) 
        sol = solve(prob)
    else 
        sol = solve(prob, saveat = dt)
    end

    t = sol.t 
    x = sol.u 

    return t, x  
end

export propagate_2Body

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

## ============================================ ##

export cart2kep
function cart2kep(cart, mu, tol=1e-20)
    # get position and velocity components
    r = cart[1:3]
    v = cart[4:6]

    # orbit energy calc
    h = cross(r, v) # specific angular momentum
    magR = norm(r)
    magV = norm(v)
    energy = 0.5 * (magV ^ 2) - (mu / magR) # vis viva equation

    # eccentricity vector
    ecc = (cross(v, h) - mu * (r / magR)) / mu
    magEcc = norm(ecc) 

    # semi-latus rectum 
    p = (norm(h) ^ 2) / mu 

    # semi-major axis 
    sma = p / (1 - magEcc^2) 

    # # determine orbit type based on eccentricity magnitude
    # if (magEcc <= 1.0)
    #     # ellipse
    #     sma = -mu / (2.0 * energy)
    # else
    #     # hyperbola
    #     sma = mu / (2.0 * energy)
    # end

    # get mean anomaly
    theta = acos(dot(r, ecc) / (magR * magEcc))
    if (dot(r, v) < 0.0)
        theta = 2.0 * np.pi - theta
    end
    inc = acos(h[3] / norm(h))

    K = [0.0; 0.0; 1.0]
    n = cross(K, h)
    normN = norm(n)

    # get right ascension
    raan = acos(n[1] / normN)
    if n[2] < 0.0
        raan = 2.0 * pi - raan
    end

    # get argument of periapsis
    omega = acos(dot(n, ecc) / (normN * magEcc))
    if ecc[3] < 0.0
        omega = 2.0 * pi - omega
    end

    # singularity checks
    I = [1.0; 0.0; 0.0]

    if magEcc < tol && inc < tol
        raan = 0.0
        omega = 0.0

        # set true longitude of periapsis as theta
        theta = acos(dot(r, I) / magR)
        if r[2] < 0.0
            theta = 2.0 * pi - theta
        end
    elseif magEcc < tol
        omega = 0.0

        # set argument of latitude as theta
        theta = acos(dot(n, r) / (normN * magR))
        if r[3] < 0.0
            theta = 2.0 * pi - theta
        end
    elseif inc < tol
        raan = 0.0

        # set longitude of periapsis as omega
        omega = acos(dot(ecc, I) / magEcc)
        if ecc[2] < 0.0
            omega = 2.0 * pi - omega
        end
    end

    # output theta here refers to true anomaly!
    kepState = [sma; magEcc; inc; raan; omega; theta]

    return kepState
end

## ============================================ ##

export orbitPeriod 
function orbitPeriod(kep, mu)
    T = 2.0 * pi * sqrt(kep[1]^3 / mu)
    return T
end