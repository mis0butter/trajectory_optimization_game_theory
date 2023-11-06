#============================================================
PROPKEPUV:

Description: Universal propagation algorithm

Inputs:
    1. x̄ₒ - Non-dimensionalized initial state vector of form [r̄; v̄]
    2. ΔUV - Universal Variable

Outputs:
    1. xf - Final state vector
    2. Δt - Change in time between initial and final states
============================================================#

function propKepUV(
    x₀_vec::AbstractVector{T}, 
    ΔUV::T
    ) where T<:Real

    # PREALLOCATING
    x0 = x₀_vec
    p = zeros(T, 7)
    
    # BREAKING UP VECTOR
    r_vec = x₀_vec[1:3]
    v_vec = x₀_vec[4:6]
    
    # SETTING CONSTANTS
    epsswitch=0.01; 
    NEGepsswitch = -epsswitch;
    onebyQfact = [1, 0.5, 0.166666666666667, 0.0416666666666667, 0.00833333333333333, 0.00138888888888889, 0.000198412698412698, 2.48015873015873e-05, 2.75573192239859e-06];
    
    # ASSIGNING UNIVERSAL VARIABLE
    x=ΔUV;
    xsq = x^2;
    
    # CALCULATING REUSED VALUES
    v0sq=v_vec⋅v_vec
    r0sq=r_vec⋅r_vec;
    r0 = sqrt(r0sq);
    onebyr0=1/r0;
    capZ = xsq * (-v0sq * r0 + 2)  * onebyr0;
    
    # CHECKING FOR OVERFLOW OF HYPERBOLICS
    if abs(capZ)>50400
        @warn "Overflow due to Excessive Hyperbolic Trajectory" capZ
        return fill(NaN, 6)
    end
    
    # HANDLING BRANCHES
    if capZ>epsswitch
        t1 = sqrt(capZ)
        t2 = cos(t1)
        t3 = 1 - t2
        t5 = -r0 * v0sq + 2
        t6 = 1 / t5
        t7 = t6 * t3
        t8 = 1 - t7
        t13 = x0[1] * x0[4] + x0[2] * x0[5] + x0[3] * x0[6]
        t15 = r0 * t6
        t16 = t15 * t3 * t13
        t18 = sin(t1)
        t22 = 0.1e1 / t1 / capZ * (t1 - t18)
        t26 = -onebyr0 * t5 * xsq * t22 + 1
        t27 = t26 * x * r0
        t28 = t16 + t27
        t42 = 0.1e1 / (t26 * x * t13 + t2 * r0 + r0 * t7)
        t43 = t42 * onebyr0
        t44 = -t26 * x
        t49 = -t15 * t3 * t42 + 1
        p[1] = x0[4] * t28 + x0[1] * t8
        p[2] = x0[5] * t28 + x0[2] * t8
        p[3] = x0[6] * t28 + x0[3] * t8
        p[4] = x0[1] * t44 * t43 + x0[4] * t49
        p[5] = x0[2] * t44 * t43 + x0[5] * t49
        p[6] = x0[3] * t44 * t43 + x0[6] * t49
        p[7] = t22 * xsq * x + t16 + t27
        if isnan(p[7])
            @warn "Something went wrong: Elliptical?" p capZ NEGepsswitch epsswitch
        end
        
    elseif capZ<NEGepsswitch
        t1 = sqrt(-capZ)
        t2 = cosh(t1)
        t3 = 1 - t2
        t5 = -r0 * v0sq + 2
        t6 = 1 / t5
        t7 = t6 * t3
        t8 = 1 - t7
        t13 = x0[1] * x0[4] + x0[2] * x0[5] + x0[3] * x0[6]
        t15 = r0 * t6
        t16 = t15 * t3 * t13
        t18 = sinh(t1)
        t22 = -0.1e1 / t1 / capZ * (t18 - t1)
        t26 = -onebyr0 * t5 * xsq * t22 + 1
        t27 = t26 * x * r0
        t28 = t16 + t27
        t42 = 0.1e1 / (t26 * x * t13 + t2 * r0 + r0 * t7)
        t43 = t42 * onebyr0
        t44 = -t26 * x
        t49 = -t15 * t3 * t42 + 1
        p[1] = x0[4] * t28 + x0[1] * t8
        p[2] = x0[5] * t28 + x0[2] * t8
        p[3] = x0[6] * t28 + x0[3] * t8
        p[4] = x0[1] * t44 * t43 + x0[4] * t49
        p[5] = x0[2] * t44 * t43 + x0[5] * t49
        p[6] = x0[3] * t44 * t43 + x0[6] * t49
        p[7] = t22 * xsq * x + t16 + t27
        if isnan(p[7])
            @warn "Something went wrong: Hyperbolic?" p capZ NEGepsswitch epsswitch
        end
        
    else
        t3 = -r0 * v0sq + 2
        t4 = t3 * xsq
        t7 = xsq ^ 2
        t8 = t3 ^ 2
        t9 = t8 * t7
        t10 = 0.1e1 / r0sq
        t13 = x ^ 2
        t14 = t13 ^ 2
        t17 = t8 * t3 * t14 * t13
        t18 = t10 * onebyr0
        t21 = 0.1e1 / 0.2e1 - onebyr0 * onebyQfact[4] * t4 + onebyQfact[6] * t10 * t9 - onebyQfact[8] * t18 * t17
        t23 = -t21 * onebyr0 * xsq + 1
        t28 = x0[1] * x0[4] + x0[2] * x0[5] + x0[3] * x0[6]
        t30 = t21 * xsq * t28
        t38 = 0.1e1 / 0.6e1 - onebyr0 * onebyQfact[5] * t4 + onebyQfact[7] * t10 * t9 - onebyQfact[9] * t18 * t17
        t40 = onebyr0 * t3
        t42 = -t40 * xsq * t38 + 1
        t43 = t42 * x * r0
        t44 = t30 + t43
        t53 = t21 * xsq
        t60 = 0.1e1 / (t53 + t42 * x * t28 + (-t40 * t53 + 1) * r0)
        t61 = t60 * onebyr0
        t62 = -t42 * x
        t67 = -t21 * t60 * xsq + 1
        p[1] = x0[1] * t23 + x0[4] * t44
        p[2] = x0[2] * t23 + x0[5] * t44
        p[3] = x0[3] * t23 + x0[6] * t44
        p[4] = x0[1] * t62 * t61 + x0[4] * t67
        p[5] = x0[2] * t62 * t61 + x0[5] * t67
        p[6] = x0[3] * t62 * t61 + x0[6] * t67
        p[7] = t38 * x * xsq + t30 + t43
        if isnan(p[7])
            @warn "Something went wrong: Parabolic?" p capZ NEGepsswitch epsswitch
            @warn "Other Values" xsq v0sq r0  onebyr0
        end
    end
    
    # RETURNING
    xf = p[1:6]
    Δt = p[7]
    return xf, Δt
    
end