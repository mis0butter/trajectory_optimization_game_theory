struct SpacecraftDynamics{T1,T2} <: AbstractDynamics
    dt::Float64
    mu::Float64
    state_bounds::T1
    control_bounds::T2
    integration_scheme::Symbol

    function SpacecraftDynamics(;
        dt = 1.0,
        mu = 398600.4415,
        state_bounds::T1 = (; lb = [-Inf, -Inf, -Inf, -Inf, -Inf, -Inf], ub = [Inf, Inf, Inf, Inf, Inf, Inf]),
        control_bounds::T2 = (; lb = [-Inf, -Inf, -Inf], ub = [Inf, Inf, Inf]),
        integration_scheme = :twoBody,
    ) where {T1,T2}
        supported_integration_schemes = (:twoBody, :withJ2, :withJ2andDrag)
        integration_scheme ∈ supported_integration_schemes ||
            throw(ArgumentError("integration_scheme must be one of $supported_integration_schemes"))
        new{T1,T2}(dt, mu, state_bounds, control_bounds, integration_scheme)
    end
end

function TrajectoryGamesBase.horizon(sys::SpacecraftDynamics)
    ∞
end

function TrajectoryGamesBase.state_dim(dynamics::SpacecraftDynamics)
    6
end

function TrajectoryGamesBase.control_dim(dynamics::SpacecraftDynamics)
    3
end

function TrajectoryGamesBase.state_bounds(dynamics::SpacecraftDynamics)
    dynamics.state_bounds
end

function TrajectoryGamesBase.control_bounds(dynamics::SpacecraftDynamics)
    dynamics.control_bounds
end

function (sys::SpacecraftDynamics)(state, control, t)
    x, y, z, vx, vy, vz = state
    dv_x, dv_y, dv_z = control
    dt = sys.dt
    mu = sys.mu
    h = dt

    # apply control input to velocity components
    vx += dv_x
    vy += dv_y
    vz += dv_z

    stateWithControl = [x, y, z, vx, vy, vz]

    # propagate with 2-body dynamics, in Kelperian 
    xKOE = cart2kep(stateWithControl, mu)
    xKOE[6] += (2*pi / orbitPeriod(xKOE, mu)) * dt

    # convert orbit back to Cartesian
    stateOut = kep2cart(xKOE, mu)

    # single step RK4 integration 
    # k1 = eom_2body(stateWithControl, mu)
    # k2 = eom_2body(stateWithControl + h/2 * k1, mu)
    # k3 = eom_2body(stateWithControl + h/2 * k2, mu)
    # k4 = eom_2body(stateWithControl + h * k3, mu)
    # state_out = stateWithControl + h/6 * (k1 + 2*k2 + 2*k3 + k4)

    x_t, y_t, z_t, vx_t, vy_t, vz_t = state_out

    # next state
    [x_t, y_t, z_t, vx_t, vy_t, vz_t]
end

# function eom_2body(state, mu)

#     f = zeros(6)
    
#     f[1] = state[4]
#     f[2] = state[5]
#     f[3] = state[6]

#     mu_div_r3 = -mu/sqrt(state[1]^2 + state[2]^2 + state[3]^2) ^ 3
#     f[4] = mu_div_r3 * state[1]
#     f[5] = mu_div_r3 * state[2]
#     f[6] = mu_div_r3 * state[3]

#     return f
# end
