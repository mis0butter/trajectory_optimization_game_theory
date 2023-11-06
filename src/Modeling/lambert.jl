#============================================================
LAMBERT:

Description: Calculates Lambert transfer between two points in space. `xᵢ_vec` is the initial state vector, 
			 `xⱼ_vec` is the final state vector.

Inputs:
	1. xᵢ_vec - initial state vector
	2. xⱼ_vec - final state vector
	3. tof - final time
	4. μ - gravitational parameter
	5. iterMax - maximum amount of iterations, initialized to 100 by default
	6. errorTol - error tolerance, initialized to 1.0e-5 by default

Outputs:
	1. Lambert structure containing the transfer
============================================================#
function lambert(
    xᵢ_vec::AbstractVector{T}, 
    xⱼ_vec::AbstractVector{T}, 
    tof::T, 
    μ::T;
    iterMax::Int = 100,
    errorTol::Float64 = 1e-5
    ) where T

	#=====================================================================
	C2C3:

	Description: Sub-routine that converts ψ to c2 and c3

	Inputs:
		1. ψ - constant 
	
	Outputs:
		1. c2 - constant
		2. c3 - constant
	=====================================================================#
	function c2c3(ψ)
		if ψ > 1e-6
			c2 = (1.0 - cos(√ψ))/ψ
			c3 = (√ψ - sin(√ψ))/√(ψ^3)
			
		elseif ψ < -1e-6
			c2 = (1.0 - cosh(√(-ψ)))/ψ
			c3 = (sinh(√(-ψ)) - √(-ψ))/√(-ψ^3)

		else
			c2 = 1/2
			c3 = 1/6

		end

		return c2, c3
	end
	
	# HANDLING INPUTS
	rᵢ_vec = xᵢ_vec[1:3]
	vᵢ_vec = xᵢ_vec[4:6]
	rᵢ = norm(rᵢ_vec)
	rⱼ_vec = xⱼ_vec[1:3]
	vⱼ_vec = xⱼ_vec[4:6]
	rⱼ = norm(rⱼ_vec)

	# FINDING POSITIONAL ANGLES
	νᵢ = atand(rᵢ_vec[2], rᵢ_vec[1])
	νⱼ = atand(rⱼ_vec[2], rⱼ_vec[1])
	δν = νⱼ - νᵢ
	δν += 360.0
	δν = mod(δν, 360.0)
	dir = δν > 180 ? -1 : 1

	# ETC
	cδν = rᵢ_vec⋅rⱼ_vec / (rᵢ*rⱼ)
	A = dir*√(rᵢ*rⱼ*(1 + cδν))

	# DEFINING INITIAL CONSTANTS
	c2 = 1/2
	c3 = 1/6
	ψ = 0.
	ψᵤ = 4*π^2
	ψₗ = -4*π

	y = rᵢ + rⱼ + ((A*(ψ*c3 - 1))/√(c2))
	Xᵢ = √(y/c2)
	δtᵢ = (Xᵢ^3 * c3 + A*√(y))/√(μ)

	N = 0.8
	iter = 0

	# CONVERGING SOLUTIONS
	while abs(δtᵢ - tof) > errorTol
		# RECALCULATING CONSTANTS
		y = rᵢ + rⱼ + ((A*(ψ*c3 - 1))/√(c2))
		if A > 0.0 && y < 0.0
			ψ = N*(1/c3)*(1 - ((√(c2)/A)*(rᵢ + rⱼ)))
			c2, c3 = c2c3(ψ)
			y = rᵢ + rⱼ + ((A*(ψ*c3 - 1))/√(c2))
		end
		if y/c2 < 0
			iter = iterMax
			break
		end
		
		Xᵢ = √(y/c2)
		δtᵢ = (Xᵢ^3 * c3 + A*√(y))/√(μ)

		if δtᵢ ≤ tof
			ψₗ = ψ
		else
			ψᵤ = ψ
		end

		ψ = 0.5*(ψₗ+ψᵤ)
		c2, c3 = c2c3(ψ)
		
		# ITERATING
		iter += 1
		if iter == iterMax; break; end
	end

	
	if iter < iterMax
		# f & g FUNCTIONS
		f = 1 - y/rᵢ
		g = A*√(y/μ)
		ġ = 1 - y/rⱼ

		# FINDING FINAL VALUES
		# vᵢ_vec = rdiv!(rⱼ_vec - lmul!(f, rᵢ_vec), g)
		vᵢ_vec = (rⱼ_vec - f*rᵢ_vec)/g
		# vⱼ_vec = rdiv!(lmul!(ġ, r_vecⱼ) - r_vecᵢ, g)
		vⱼ_vec = (ġ*rⱼ_vec - rᵢ_vec)/g
		Δvᵢ = xᵢ_vec[4:6]-vᵢ_vec
		Δvⱼ = xⱼ_vec[4:6]-vⱼ_vec

	else
		Δvᵢ = [10., 10, 10]
		Δvⱼ = [10., 10, 10]
	end

	# OUTPUTTING
	return Lambert(
        vcat(rᵢ_vec, vᵢ_vec), 
        vcat(rⱼ_vec, vⱼ_vec), 
        δtᵢ, 
        Δvᵢ, 
        Δvⱼ, 
        iter,
        δtᵢ - tof,
        0
    )
end