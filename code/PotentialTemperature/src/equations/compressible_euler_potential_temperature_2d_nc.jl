using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, FluxLMARS, boundary_condition_slip_wall, have_nonconservative_terms

@muladd begin
#! format: noindent
struct CompressibleEulerPotentialTemperatureEquations2DNC{RealT <: Real} <:
	   AbstractCompressibleEulerEquations{2, 5}
	p_0::RealT
	c_p::RealT
	c_v::RealT
	g::RealT
	R::RealT
	gamma::RealT
	inv_gamma_minus_one::RealT
end

function CompressibleEulerPotentialTemperatureEquations2DNC(; g = 9.81, RealT = Float64)
	p_0 = 100_000.0
	c_p = 1004.0
	c_v = 717.0
	R = c_p - c_v
	gamma = c_p / c_v
	inv_gamma_minus_one = inv(gamma - 1.0)
	return CompressibleEulerPotentialTemperatureEquations2DNC{RealT}(p_0, c_p, c_v, g, R,
		gamma, inv_gamma_minus_one)
end

function varnames(::typeof(cons2cons),
	::CompressibleEulerPotentialTemperatureEquations2DNC)
	("rho", "rho_v1", "rho_v2", "rho_theta", "phi")
end

varnames(::typeof(cons2prim), ::CompressibleEulerPotentialTemperatureEquations2DNC) = ("rho",
	"v1",
	"v2",
	"p1", "phi")

have_nonconservative_terms(::CompressibleEulerPotentialTemperatureEquations2DNC) = Trixi.True()


# Calculate 1D flux for a single point in the normal direction.
# Note, this directional vector is not normalized.
@inline function flux(u, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, phi = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
	rho_v_normal = rho * v_normal
	f1 = rho_v_normal
	f2 = (rho_v_normal) * v1 + p * normal_direction[1]
	f3 = (rho_v_normal) * v2 + p * normal_direction[2]
	f4 = (rho_theta) * v_normal
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux(u, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, phi = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma

	if orientation == 1
		f1 = rho_v1
		f2 = rho_v1 * v1 + p
		f3 = rho_v1 * v2
		f4 = rho_theta * v1
	else
		f1 = rho_v2
		f2 = rho_v2 * v1
		f3 = rho_v2 * v2 + p
		f4 = rho_theta * v2
	end

	return SVector(f1, f2, f3, f4, 0)
end

# Slip-wall boundary condition
# Determine the boundary numerical surface flux for a slip wall condition.
# Imposes a zero normal velocity at the wall.
@inline function boundary_condition_slip_wall_2(u_inner, normal_direction::AbstractVector,
	x, t,
	surface_flux_functions,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# normalize the outward pointing direction
	normal = normal_direction / norm(normal_direction)
	surface_flux_function, nonconservative_flux_function = surface_flux_functions

	# compute the normal velocity
	u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

	# create the "external" boundary solution state
	u_boundary = SVector(u_inner[1],
		u_inner[2] - 2 * u_normal * normal[1],
		u_inner[3] - 2 * u_normal * normal[2],
		u_inner[4], u_inner[5])

	# calculate the boundary flux
	flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
	noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
		equations)
	return flux, noncons_flux
end

@inline function boundary_condition_slip_wall_2(u_inner, orientation,
	direction, x, t,
	surface_flux_functions,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	surface_flux_function, nonconservative_flux_function = surface_flux_functions

	## get the appropriate normal vector from the orientation
	if orientation == 1
		u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4], u_inner[5])
	else # orientation == 2
		u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4], u_inner[5])
	end
	surface_flux_function = surface_flux_functions[1]
	# Calculate boundary flux
	if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
		flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
		noncons_flux = nonconservative_flux_function(u_inner, u_boundary, orientation,
			equations)
	else # u_boundary is "left" of boundary, u_inner is "right" of boundary
		flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
		noncons_flux = nonconservative_flux_function(u_boundary, u_inner, orientation,
			equations)
	end

	return flux, noncons_flux
end

@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
	x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	@unpack gamma = equations
	norm_ = norm(normal_direction)
	# Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
	normal = normal_direction / norm_

	# rotate the internal solution state
	u_local = rotate_to_x(u_inner, normal, equations)

	# compute the primitive variables
	rho_local, v_normal, v_tangent, p_local, phi = cons2prim(u_local, equations)
	# Get the solution of the pressure Riemann problem
	# See Section 6.3.3 of
	# Eleuterio F. Toro (2009)
	# Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
	# [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
	if v_normal <= 0.0
		sound_speed = sqrt(gamma * p_local / rho_local) # local sound speed
		p_star = p_local *
				 (1.0 + 0.5f0 * (gamma - 1) * v_normal / sound_speed)^(2.0 * gamma * equations.inv_gamma_minus_one)
	else # v_normal > 0.0
		A = 2.0 / ((gamma + 1) * rho_local)
		B = p_local * (gamma - 1) / (gamma + 1)
		p_star = p_local +
				 0.5f0 * v_normal / A *
				 (v_normal + sqrt(v_normal^2 + 4.0 * A * (p_local + B)))
	end

	# For the slip wall we directly set the flux as the normal velocity is zero
	return SVector(zero(eltype(u_inner)),
		p_star * normal[1],
		p_star * normal[2],
		zero(eltype(u_inner)), zero(eltype(u_inner))) * norm_
end

# Fix sign for structured mesh.
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
	direction, x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
	# to be inward pointing on the -x and -y sides due to the orientation convention used by StructuredMesh
	if isodd(direction)
		boundary_flux = -boundary_condition_slip_wall(u_inner, -normal_direction,
			x, t, surface_flux_function,
			equations)
	else
		boundary_flux = boundary_condition_slip_wall(u_inner, normal_direction,
			x, t, surface_flux_function,
			equations)
	end

	return boundary_flux
end

@inline function flux_nonconservative_gravity_am(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, _, _ = u_ll
	rho_rr, _, _, _, _ = u_rr
	rho_avg = (rho_ll + rho_rr) * 0.5f0

	phi_jump = u_rr[5] - u_ll[5]
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	return SVector(0,
		0,
		normal_direction[2] * 9.81 * rho_avg * phi_jump
		#normal_direction[2]*9.81*rho_ll*u_rr[5]
		, 0, 0)
end

@inline function flux_nonconservative_gravity_log(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, _, _ = u_ll
	rho_rr, _, _, _, _ = u_rr
	rho_avg = ln_mean(rho_ll, rho_rr)
	RealT = eltype(u_ll)
	phi_jump = u_rr[5] - u_ll[5]
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	return SVector(0,
		normal_direction[1] * RealT(9.81) * rho_avg * phi_jump,
		normal_direction[2] * RealT(9.81) * rho_avg * phi_jump
		#normal_direction[2]*9.81*rho_ll*u_rr[5]
		, 0, 0)
end

@inline function flux_nonconservative_gravity_gamma(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, _, _ = u_ll
	rho_rr, _, _, _, _ = u_rr

	phi_jump = u_rr[5] - u_ll[5]
	rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
	#rho_avg = ln_mean(rho_ll, rho_rr)
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	if orientation == 1
		f = SVector(0, 9.81 * rho_avg * phi_jump, 0, 0, 0)
	else # orientation == 2
		f = SVector(0, 0, 9.81 * rho_avg * phi_jump, 0, 0)
		#f = SVector(0, 0, 9.81 * rho_avg * u_rr[5], 0, 0)
	end
	return f
end

@inline function flux_nonconservative_gravity_gamma(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, _, _ = u_ll
	rho_rr, _, _, _, _ = u_rr
	rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
	#rho_avg = ln_mean(rho_ll, rho_rr)
	phi_jump = u_rr[5] - u_ll[5]
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	return SVector(0,
		normal_direction[1] * 9.81 * rho_avg * phi_jump,
		normal_direction[2] * 9.81 * rho_avg * phi_jump
		#normal_direction[2]*9.81*rho_ll*u_rr[5]
		, 0, 0)
end

@inline function initial_condition_density_wave(x, t, equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	RealT = eltype(x)
	v1 = convert(RealT, 0.1)
	v2 = convert(RealT, 0.2)
	rho = 1 + convert(RealT, 0.98) * sinpi(2 * (x[1] + x[2] - t * (v1 + v2)))
	rho_v1 = rho * v1
	rho_v2 = rho * v2
	p = 20
	rho_e = p / (equations.gamma - 1) + 0.5f0 * rho * (v1^2 + v2^2)
	return prim2cons(SVector(rho, v1, v2, p, 0), equations)
end

# Rotate momentum flux. The same as in compressible Euler.
@inline function rotate_to_x(u, normal_vector::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# cos and sin of the angle between the x-axis and the normalized normal_vector are
	# the normalized vector's x and y coordinates respectively (see unit circle).
	c = normal_vector[1]
	s = normal_vector[2]

	# Apply the 2D rotation matrix with normal and tangent directions of the form
	# [ 1    0    0   0;
	#   0   n_1  n_2  0;
	#   0   t_1  t_2  0;
	#   0    0    0   1 ]
	# where t_1 = -n_2 and t_2 = n_1

	return SVector(u[1],
		c * u[2] + s * u[3],
		-s * u[2] + c * u[3],
		u[4], u[5])
end

# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	a = flux_lmars.speed_of_sound
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	norm_ = norm(normal_direction)

	rho = 0.5f0 * (rho_ll + rho_rr)

	p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_rr - v_ll) / norm_
	v_interface = 0.5f0 * (v_ll + v_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

	if (v_interface > 0)
		f1, f2, f3, f4 = u_ll * v_interface
	else
		f1, f2, f3, f4 = u_rr * v_interface
	end

	return SVector(f1,
		f2 + p_interface * normal_direction[1],
		f3 + p_interface * normal_direction[2],
		f4, 0)
end

## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation
@inline function flux_theta_surface(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	rho_rr_wb = rho_rr^equations.gamma
	rho_ll_wb = rho_ll^equations.gamma
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	K = equations.p_0 * (equations.R / equations.p_0)^equations.gamma
	#if normal_direction[1] == 0
	p_avg = 0.5f0 * (p_ll + p_rr)
	#p_avg2 = 0.5f0 * (rho_rr_wb + rho_ll_wb) * K
	#@show p_avg
	#p_avg = 0.5f0 * ((rho_ll*300)^equations.gamma + (rho_rr*300)^equations.gamma) * K
	#p_avg = 0.5f0 * 300^equations.gamma*((rho_ll)^equations.gamma + rho_rr^equations.gamma)*K
	#@show p_avg
	#@show normal_direction
	#p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	#@show f1, f2, f4
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_theta(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	rho_rr_wb = rho_rr^equations.gamma
	rho_ll_wb = rho_ll^equations.gamma
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	K = equations.p_0 * (equations.R / equations.p_0)^equations.gamma
	#if normal_direction[1] == 0
	p_avg = 0.5f0 * (p_ll + p_rr)
	#p_avg2 = 0.5f0 * (rho_rr_wb + rho_ll_wb) * K
	#@show p_avg
	#p_avg = 0.5f0 * ((rho_ll*300)^equations.gamma + (rho_rr*300)^equations.gamma) * K
	#p_avg = 0.5f0 * 300^equations.gamma*((rho_ll)^equations.gamma + rho_rr^equations.gamma)*K
	#@show p_avg
	#@show normal_direction
	#p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	#@show f1, f2, f4
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_theta_mod(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	rho_mean = 0.5 * (rho_ll + rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = 0.50f0 * (rho_theta_ll + rho_theta_rr) * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_theta_sto(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)

	rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	rho_theta_ll = last(u_ll)
	rho_theta_rr = last(u_rr)

	# Compute the necessary mean values
	rho_mean = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_theta_global(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	#rho_mean = ln_mean(rho_ll, rho_rr)
	#rho_mean = 0.5f0 * (rho_ll + rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	#f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f1 = f4 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_theta_rhos(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	return SVector(f1, f2, f3, f4, 0)
end

@inline function boundary_condition_slip_wall2(u_inner, normal_direction::AbstractVector,
	x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# normalize the outward pointing direction
	normal = normal_direction / norm(normal_direction)

	# compute the normal velocity
	u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

	# create the "external" boundary solution state
	u_boundary = SVector(u_inner[1],
		u_inner[2] - 2 * u_normal * normal[1],
		u_inner[3] - 2 * u_normal * normal[2],
		u_inner[4], u_inner[5])

	# calculate the boundary flux
	flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

	return flux
end

"""
boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
surface_flux_function, equations::ShallowWaterEquations2D)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall2(u_inner, orientation,
	direction, x, t,
	surface_flux_functions,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	surface_flux_function, nonconservative_flux_function = surface_flux_functions
	## get the appropriate normal vector from the orientation
	if orientation == 1
		u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4], u_inner[5])
	else # orientation == 2
		u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4], u_inner[5])
	end

	# Calculate boundary flux
	if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
		flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
	else # u_boundary is "left" of boundary, u_inner is "right" of boundary
		flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
	end

	return flux
end

@inline function flux_theta_AM(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5 * (rho_ll + rho_rr)

	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_theta(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	println("Testing time")
	@time gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)
	@time gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = gammamean * v1_avg
	else
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = gammamean * v2_avg
	end

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_theta_AM(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5 * (rho_ll + rho_rr)

	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = gammamean * v1_avg
	else
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = gammamean * v2_avg
	end

	return SVector(f1, f2, f3, f4, 0)
end

# Entropy stable, density and pressure positivity preserving flux
@inline function flux_theta_es(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	_, v1_ll, v2_ll, _, _ = cons2prim(u_ll, equations)
	_, v1_rr, v2_rr, _, _ = cons2prim(u_rr, equations)

	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = gammamean * v1_avg
	else
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = gammamean * v2_avg
	end
	if orientation == 1
		lambda = max(abs(v1_ll), abs(v1_rr))
	else
		lambda = max(abs(v2_ll), abs(v2_rr))
	end

	f1 = f1 - lambda * 0.5 * (rho_rr - rho_ll)
	f2 = f2 - lambda * 0.5 * (rho_rr * v1_rr - rho_ll * v1_ll)
	f3 = f3 - lambda * 0.5 * (rho_rr * v2_rr - rho_ll * v2_ll)
	f4 = f4 - lambda * 0.5 * (rho_theta_rr - rho_theta_ll)
	return SVector(f1, f2, f3, f4, 0)
end


@inline function flux_theta_es(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	_, _, _, rho_theta_ll, _ = u_ll
	_, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	lambda = max(abs(v_dot_n_ll), abs(v_dot_n_rr))


	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) - 0.5f0 * lambda * u_rr[1] - u_ll[1]
	f2 = f1 * v1_avg + p_avg * normal_direction[1] - 0.5f0 * lambda * u_rr[2] - u_ll[2]
	f3 = f1 * v2_avg + p_avg * normal_direction[2] - 0.5f0 * lambda * u_rr[3] - u_ll[3]
	f4 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) - 0.5f0 * lambda * u_rr[4] - u_ll[4]
	return SVector(f1, f2, f3, f4, 0)
end

@inline function prim2cons(prim,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, v1, v2, p, phi = prim
	rho_v1 = rho * v1
	rho_v2 = rho * v2
	rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
	return SVector(rho, rho_v1, rho_v2, rho_theta, phi)
end

@inline function cons2prim(u,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, phi = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	return SVector(rho, v1, v2, p, phi)
end

@inline function cons2cons(u,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	return u
end

@inline function cons2entropy(u,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, _ = u

	k = equations.p_0 * (equations.R / equations.p_0)^equations.gamma
	w1 = -0.5f0 * rho_v1^2 / (rho)^2 - 0.5f0 * rho_v2^2 / (rho)^2
	w2 = rho_v1 / rho
	w3 = rho_v2 / rho
	w4 = equations.gamma / (equations.gamma - 1) * k * (rho_theta)^(equations.gamma - 1)

	return SVector(w1, w2, w3, w4, 0)
end

@inline function entropy_math(cons,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	# Mathematical entropy
	-

	U = (p / (equations.gamma - 1) + 0.5f0 * (cons[2]^2 + cons[3]^2) / (cons[1]))

	return U
end

# Default entropy is the mathematical entropy
@inline function entropy(cons,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	entropy_math(cons, equations)
end

@inline function energy_total(cons,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	entropy(cons, equations)
end

@inline function energy_kinetic(cons,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	return 0.5f0 * (cons[2]^2 + cons[3]^2) / (cons[1])
end

@inline function max_abs_speeds(u,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, v1, v2, p, _ = cons2prim(u, equations)
	c = sqrt(equations.gamma * p / rho)

	return abs(v1) + c, abs(v2) + c
end

@inline function density_pressure(u, equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, _ = u
	rho_times_p = rho * equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	return rho_times_p
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	# Get the velocity value in the appropriate direction
	if orientation == 1
		v_ll = v1_ll
		v_rr = v1_rr
	else # orientation == 2
		v_ll = v2_ll
		v_rr = v2_rr
	end
	# Calculate sound speeds
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	λ_max = max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function stolarsky_mean_opt(x::RealT, y::RealT, gamma::RealT, xg::RealT, yg::RealT) where {RealT <: Real}
	epsilon_f2 = convert(RealT, 1.0e-4)
	f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
	if f2 < epsilon_f2
		# convenience coefficients
		c1 = convert(RealT, 1 / 3) * (gamma - 2)
		c2 = convert(RealT, -1 / 15) * (gamma + 1) * (gamma - 3) * c1
		c3 = convert(RealT, -1 / 21) * (2 * gamma * (gamma - 2) - 9) * c2
		return 0.5f0 * (x + y) * @evalpoly(f2, 1, c1, c2, c3)
	else
		return (gamma - 1) / gamma * (yg - xg) /
			   (yg / y - xg / x)
	end
end

@inline function well_balanced_v1(u, equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, _ = u
	v1 = rho_v1 / rho
	return abs(v1)
end

@inline function well_balanced_v2(u, equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, _ = u
	v2 = rho_v2 / rho
	return abs(v2)
end

end # @muladd
