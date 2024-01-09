using trajectory_optimization_game_theory
using LinearAlgebra 

## ============================================ ## 
# define IC, target state, and lambert solve 

r_0     = [20.0e6, 20.0e6, 0]   # [m] 
r_f     = [-20.0e6, 10.0e6, 0]  # [m] 
tof     = 1.0 * 86400 
mu      = 398600.4418e9         # [m^3/s^2] 
dm      = "pro" 
Dtsec   = tof 

## ============================================ ##
# solve and propagate lambert orbit 

v_0, v_f = lambertbattin( r_0, r_f, mu, dm, tof ) 

rv_0 = [ r_0; v_0 ] 

## ============================================ ##
# cart2kep 

cart = [ r_0 ; v_0*1.3 ]
kep  = cart2kep( cart, mu ) 

# cart = rv_0 

    # get position and velocity components
    r = cart[1:3]
    v = cart[4:6]

    # orbit energy calc
    h = cross(r, v) # specific angular momentum
    r_mag = norm(r)
    v_mag = norm(v)
    energy = 0.5 * (v_mag ^ 2) - (mu / r_mag) # vis viva equation

    # eccentricity vector
    ecc = (cross(v, h) - mu * (r / r_mag)) / mu
    magEcc = norm(ecc)

    # determine orbit type based on eccentricity magnitude
    if (magEcc <= 1.0)
        # ellipse
        sma = -mu / (2.0 * energy)
    else
        # hyperbola
        sma = mu / (2.0 * energy)
    end

    # get mean anomaly
    theta = acos(dot(r, ecc) / (r_mag * magEcc))
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
        theta = acos(dot(r, I) / r_mag)
        if r[2] < 0.0
            theta = 2.0 * pi - theta
        end
    elseif magEcc < tol
        omega = 0.0

        # set argument of latitude as theta
        theta = acos(dot(n, r) / (normN * r_mag))
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

## ============================================ ##

rv = [ r_0 ; v_0*1.3 ]

# my function 
# function rv2oe( rv, mu = 1.0 ) 

    # get position and velocity components
    r_vec = rv[1:3]
    v_vec = rv[4:6]

    # orbit energy calc
    h_vec = cross(r_vec, v_vec) # specific angular momentum
    r = norm(r_vec)
    v = norm(v_vec)

    # node vector 
    K_vec = [0.0; 0.0; 1.0] 
    n_vec = cross( K_vec, h_vec ) 

    # eccentricity vector 
    e_vec = 1/mu * ( ( v^2 - mu/r )*r_vec - dot(r_vec, v_vec)*v_vec )

    # semi-latus rectum 
    p = norm(h_vec)^2 / mu 

    # semi-major axis 
    a = p / ( 1 - norm(e_vec)^2 ) 

    # eccentricity 
    e = norm(e_vec) 

    # inclination 
    i = acos( h_vec[3] / norm(h_vec) ) 

    # right ascension of the ascending node 
    Omega = acos( n_vec[1] / norm(n_vec) ) 
    if n_vec[2] < 0.0 
        Omega = 2*pi - Omega 
    end 

    # argument of periapsis 
    w = acos( dot(n_vec, e_vec) / ( norm(n_vec)*norm(e_vec) ) ) 
    if e_vec[3] < 0.0 
        w = 2*pi - w 
    end 

    # true anomaly 
    nu = acos( dot(e_vec, r_vec) / ( norm(e_vec)*norm(r_vec) ) ) 
    