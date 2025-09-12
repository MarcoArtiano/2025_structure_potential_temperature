using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, FluxLMARS, boundary_condition_slip_wall, have_nonconservative_terms

@muladd begin
#! format: noindent
struct CompressibleEulerPotentialTemperatureEquations3DNC{RealT <: Real} <:
	   AbstractCompressibleEulerEquations{3, 6}
	p_0::RealT
	c_p::RealT
	c_v::RealT
	g::RealT
	R::RealT
	gamma::RealT
	inv_gamma_minus_one::RealT
	K::RealT
	stolarsky_factor::RealT
end

function CompressibleEulerPotentialTemperatureEquations3DNC(; g = 9.81, RealT = Float64)
	p_0 = 100_000.0
	c_p = 1004.0
	c_v = 717.0
	R = c_p - c_v
	gamma = c_p / c_v
	inv_gamma_minus_one = inv(gamma - 1)
	K = p_0 * (R / p_0)^gamma
	stolarsky_factor = (gamma - 1.0) / gamma
	return CompressibleEulerPotentialTemperatureEquations3DNC{RealT}(p_0, c_p, c_v, g, R,
		gamma, inv_gamma_minus_one, K, stolarsky_factor)
end

function varnames(::typeof(cons2cons),
	::CompressibleEulerPotentialTemperatureEquations3DNC)
	("rho", "rho_v1", "rho_v2", "rho_v3", "rho_theta", "phi")
end
Trixi.have_nonconservative_terms(::CompressibleEulerPotentialTemperatureEquations3DNC) = Trixi.True()

varnames(::typeof(cons2prim), ::CompressibleEulerPotentialTemperatureEquations3DNC) = ("rho",
	"v1",
	"v2",
	"v3",
	"p1", "phi")

# Calculate 1D flux for a single point in the normal direction.
# Note, this directional vector is not normalized.
@inline function flux(u, normal_direction::AbstractVector, equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, rho_v1, rho_v2, rho_v3, rho_theta, _ = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v3 = rho_v3 / rho
	rho, _, _, _, p, _ = cons2prim(u, equations)
	v_normal = v1 * normal_direction[1] + v2 * normal_direction[2] + v3 * normal_direction[3]
	rho_v_normal = rho * v_normal
	f1 = rho_v_normal
	f2 = rho_v_normal * v1 + p * normal_direction[1]
	f3 = rho_v_normal * v2 + p * normal_direction[2]
	f4 = rho_v_normal * v3 + p * normal_direction[3]
	f5 = rho_theta * v_normal
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function flux(u, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, rho_v1, rho_v2, rho_v3, rho_theta, _ = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v3 = rho_v3 / rho
	p = equations.K * rho_theta^equations.gamma
	if orientation == 1
		f1 = rho_v1
		f2 = rho_v1 * v1 + p
		f3 = rho_v1 * v2
		f4 = rho_v1 * v3
		f5 = rho_theta * v1
	elseif orientation == 2
		f1 = rho_v2
		f2 = rho_v2 * v1
		f3 = rho_v2 * v2 + p
		f4 = rho_v2 * v3
		f5 = rho_theta * v2
	else
		f1 = rho_v3
		f2 = rho_v3 * v1
		f3 = rho_v3 * v2
		f4 = rho_v3 * v3 + p
		f5 = rho_theta * v3
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
	x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	norm_ = Trixi.norm(normal_direction)
	# Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
	normal = normal_direction / norm_

	# Some vector that can't be identical to normal_vector (unless normal_vector == 0)
	tangent1 = SVector(normal_direction[2], normal_direction[3], -normal_direction[1])
	# Orthogonal projection
	tangent1 -= Trixi.dot(normal, tangent1) * normal
	tangent1 = Trixi.normalize(tangent1)

	# Third orthogonal vector
	tangent2 = Trixi.normalize(Trixi.cross(normal_direction, tangent1))

	# rotate the internal solution state
	u_local = rotate_to_x(u_inner, normal, tangent1, tangent2, equations)

	# compute the primitive variables
	rho_local, v_normal, v_tangent1, v_tangent2, p_local, _ = cons2prim(u_local, equations)

	# Get the solution of the pressure Riemann problem
	# See Section 6.3.3 of
	# Eleuterio F. Toro (2009)
	# Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction
	# [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
	if v_normal <= 0.0
		sound_speed = sqrt(equations.gamma * p_local / rho_local) # local sound speed
		p_star = p_local *
				 (1 + 0.5 * (equations.gamma - 1) * v_normal / sound_speed)^(2 *
																			 equations.gamma *
																			 equations.inv_gamma_minus_one)
	else # v_normal > 0.0
		A = 2 / ((equations.gamma + 1) * rho_local)
		B = p_local * (equations.gamma - 1) / (equations.gamma + 1)
		p_star = p_local +
				 0.5 * v_normal / A *
				 (v_normal + sqrt(v_normal^2 + 4 * A * (p_local + B)))
	end

	# For the slip wall we directly set the flux as the normal velocity is zero
	return SVector(zero(eltype(u_inner)),
		p_star * normal[1],
		p_star * normal[2],
		p_star * normal[3],
		zero(eltype(u_inner)), zero(eltype(u_inner))) * norm_
end

"""
boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
surface_flux_function, equations::CompressibleEulerEquations3D)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, orientation,
	direction, x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# get the appropriate normal vector from the orientation
	if orientation == 1
		normal_direction = SVector(1.0, 0.0, 0.0)
	elseif orientation == 2
		normal_direction = SVector(0.0, 1.0, 0.0)
	else # orientation == 3
		normal_direction = SVector(0.0, 0.0, 1.0)
	end

	# compute and return the flux using `boundary_condition_slip_wall` routine above
	return boundary_condition_slip_wall(u_inner, normal_direction, direction,
		x, t, surface_flux_function, equations)
end

"""
boundary_condition_slip_wall(u_inner, normal_direction, direction, x, t,
surface_flux_function, equations::CompressibleEulerEquations3D)

Should be used together with [`StructuredMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
	direction, x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# flip sign of normal to make it outward pointing, then flip the sign of the normal flux back
	# to be inward pointing on the -x, -y, and -z sides due to the orientation convention used by StructuredMesh
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

@inline function boundary_condition_slip_wall_2(u_inner, normal_direction::AbstractVector,
	x, t,
	surface_flux_functions,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# normalize the outward pointing direction
	normal = normal_direction / Trixi.norm(normal_direction)
	surface_flux_function, nonconservative_flux_function = surface_flux_functions
	# compute the normal velocity
	u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] + normal[3] * u_inner[4]

	# create the "external" boundary solution state
	u_boundary = SVector(u_inner[1],
		u_inner[2] - 2 * u_normal * normal[1],
		u_inner[3] - 2 * u_normal * normal[2],
		u_inner[4] - 2 * u_normal * normal[3],
		u_inner[5], u_inner[6])

	# calculate the boundary flux
	flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)
	noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
		equations)
	return flux, noncons_flux
end

function boundary_condition_slip_wall_2(u_inner, orientation, direction, x, t,
	surface_flux_functions,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Boundary state is equal to the inner state except for the velocity. For boundaries
	# in the -x/+x direction, we multiply the velocity in the x direction by -1.
	# Similarly, for boundaries in the -y/+y or -z/+z direction, we multiply the
	# velocity in the y or z direction by -1
	surface_flux_function, nonconservative_flux_function = surface_flux_functions

	if direction in (1, 2) # x direction
		u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4],
			u_inner[5], u_inner[6])
	elseif direction in (3, 4) # y direction
		u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4],
			u_inner[5], u_inner[6])
	else # z direction = (5, 6)
		u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3], -u_inner[4],
			u_inner[5], u_inner[6])
	end

	# Calculate boundary flux depending on the orientation of the boundary
	# Odd directions are in negative coordinate direction, 
	# even directions are in positive coordinate direction.
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

@inline function flux_nonconservative_gravity_am(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, _, _, _ = u_ll
	rho_rr, _, _, _, _, _ = u_rr

	phi_jump = u_rr[6] - u_ll[6]
	rho_avg = 0.5 * (rho_ll + rho_rr)
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	if orientation == 1
		f = SVector(0, 0, 0, 0, 0)
	else # orientation == 2
		f = SVector(0, 0, rho_avg * phi_jump, 0, 0)

	end
	return f
end

@inline function flux_nonconservative_gravity_am(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, _, _ = u_ll
	rho_rr, _, _, _, _ = u_rr
	rho_avg = 0.5 * (rho_ll + rho_rr)

	phi_jump = u_rr[6] - u_ll[6]
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	return SVector(0,
		normal_direction[1] * rho_avg * phi_jump,
		normal_direction[2] * rho_avg * phi_jump,
		normal_direction[3] * rho_avg * phi_jump,
		0, 0)
end

@inline function flux_nonconservative_gravity_log(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, _, _ = u_ll
	rho_rr, _, _, _, _ = u_rr
	rho_avg = ln_mean(rho_ll, rho_rr)

	phi_jump = u_rr[6] - u_ll[6]
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	return SVector(0,
		normal_direction[1] * rho_avg * phi_jump,
		normal_direction[2] * rho_avg * phi_jump,
		normal_direction[3] * rho_avg * phi_jump,
		0, 0)
end

@inline function rotate_to_x(u, normal_vector, tangent1, tangent2,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Multiply with [ 1   0        0       0   0;
	#                 0   ―  normal_vector ―   0;
	#                 0   ―    tangent1    ―   0;
	#                 0   ―    tangent2    ―   0;
	#                 0   0        0       0   1 ]
	return SVector(u[1],
		normal_vector[1] * u[2] + normal_vector[2] * u[3] +
		normal_vector[3] * u[4],
		tangent1[1] * u[2] + tangent1[2] * u[3] + tangent1[3] * u[4],
		tangent2[1] * u[2] + tangent2[2] * u[3] + tangent2[3] * u[4],
		u[5], u[6])
end


# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	a = flux_lmars.speed_of_sound
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]

	norm_ = Trixi.norm(normal_direction)

	rho = 0.5f0 * (rho_ll + rho_rr)

	p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v_rr - v_ll) / norm_
	v_interface = 0.5f0 * (v_ll + v_rr) - 1 / (2 * a * rho) * (p_rr - p_ll) * norm_

	if (v_interface > 0)
		f1, f2, f3, f4, f5, f6 = u_ll * v_interface
	else
		f1, f2, f3, f4, f5, f6 = u_rr * v_interface
	end

	return SVector(f1,
		f2 + p_interface * normal_direction[1],
		f3 + p_interface * normal_direction[2],
		f4 + p_interface * normal_direction[3],
		f5, zero(eltype(u_ll)))
end

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	c = flux_lmars.speed_of_sound

	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	if orientation == 1
		v_ll = v1_ll
		v_rr = v1_rr
	elseif orientation == 2
		v_ll = v2_ll
		v_rr = v2_rr
	else # orientation == 3
		v_ll = v3_ll
		v_rr = v3_rr
	end

	rho = 0.5 * (rho_ll + rho_rr)
	p = 0.5 * (p_ll + p_rr) - 0.5 * c * rho * (v_rr - v_ll)
	v = 0.5 * (v_ll + v_rr) - 1 / (2 * c * rho) * (p_rr - p_ll)

	# We treat the energy term analogous to the potential temperature term in the paper by
	# Chen et al., i.e. we use p_ll and p_rr, and not p
	if v >= 0
		f1, f2, f3, f4, f5, f6 = v * u_ll
	else
		f1, f2, f3, f4, f5, f6 = v * u_rr
	end

	if orientation == 1
		f2 += p
	elseif orientation == 2
		f3 += p
	else # orientation == 3
		f4 += p
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

@inline function flux_tec_am(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)


	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

@inline function flux_tec_am(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr
	theta_mean = 0.5f0 * (theta_ll + theta_rr)

	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v1_avg
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v2_avg
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = gammamean * v3_avg
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end
@inline function flux_kennedy_gruber2(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr
	theta_mean = 0.5f0 * (theta_ll + theta_rr)
	rho_theta_mean = 0.5f0 * (rho_theta_ll + rho_theta_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) * theta_mean
	f5 = rho_theta_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function flux_kennedy_gruber2(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr
	theta_mean = 0.5f0 * (theta_ll + theta_rr)

	rho_theta_mean = 0.5f0 * (rho_theta_ll + rho_theta_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = rho_theta_mean * v1_avg
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = rho_theta_mean * v2_avg
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = rho_theta_mean * v3_avg
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

@inline function flux_kennedy_gruber(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr
	theta_mean = 0.5f0 * (theta_ll + theta_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr) * theta_mean
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

@inline function flux_kennedy_gruber(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	theta_ll = rho_theta_ll / rho_ll
	theta_rr = rho_theta_rr / rho_rr
	theta_mean = 0.5f0 * (theta_ll + theta_rr)


	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = rho_mean * theta_mean * v1_avg
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = rho_mean * theta_mean * v2_avg
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = rho_mean * theta_mean * v3_avg
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation

@inline function flux_theta_p(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)
	p_avg = Trixi.stolarsky_mean(p_ll, p_rr, 3)
	#   @show "Stolarsky"
	#   @show p_avg

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	#   @show "AVG"
	#   @show p_avg

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

@inline function flux_theta_p(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)
	p_avg = Trixi.stolarsky_mean(p_ll, p_rr, 1.2)
	#    @show "Stolarsky"
	#    @show p_avg

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	#	p_avg = 0.5f0 * (p_ll + p_rr)
	#    @show "AVG"
	#    @show p_avg
	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v1_avg
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v2_avg
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = gammamean * v3_avg
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation

@inline function flux_theta(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

@inline function flux_theta(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v1_avg
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v2_avg
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = gammamean * v3_avg
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u_ll)))
end

# Entropy stable, density and pressure positivity preserving flux
@inline function flux_theta_es(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	_, v1_ll, v2_ll, v3_ll, _, _ = cons2prim(u_ll, equations)
	_, v1_rr, v2_rr, v3_rr, _, _ = cons2prim(u_rr, equations)

	f_ec = flux_theta(u_ll, u_rr, orientation, equations)
	if orientation == 1
		lambda = max(abs(v1_ll), abs(v1_rr))
	elseif orientation == 2
		lambda = max(abs(v2_ll), abs(v2_rr))
	else
		lambda = max(abs(v3_ll), abs(v3_rr))
	end
	return f_ec - 0.5f0 * lambda * (u_rr - u_ll)
end

# Entropy stable, density and pressure positivity preserving flux
@inline function flux_theta_es(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	_, v1_ll, v2_ll, v3_ll, _, _ = cons2prim(u_ll, equations)
	_, v1_rr, v2_rr, v3_rr, _, _ = cons2prim(u_rr, equations)

	f_ec = flux_theta(u_ll, u_rr, normal_direction, equations)

	lambda = max(abs(v1_ll), abs(v1_rr)) * normal_direction[1] + max(abs(v2_ll), abs(v2_rr)) * normal_direction[2] + max(abs(v3_ll), abs(v3_rr)) * normal_direction[3]

	return f_ec - 0.5f0 * lambda * (u_rr - u_ll)
end

@inline function flux_theta_AM(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (V3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function flux_theta_AM(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_rr + rho_ll)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v1_avg
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = gammamean * v2_avg
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = gammamean * v3_avg
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation

@inline function flux_theta_rhos(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = f1 * Trixi.inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function flux_theta_rhos(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)


	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = f1 * Trixi.inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = f1 * Trixi.inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = f1 * Trixi.inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation

@inline function flux_theta_rhos_AM(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (V3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)


	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]
	f5 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function flux_theta_rhos_AM(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = 0.5f0 * (rho_ll + rho_rr)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
		f5 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	elseif orientation == 2
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
		f5 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	else
		f1 = rho_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
		f5 = f1 * inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation

@inline function flux_theta_global(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2] + v3_ll * normal_direction[3]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2] + v3_rr * normal_direction[3]
	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr

	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f5 = gammamean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f1 = f5 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * v3_avg + p_avg * normal_direction[3]

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function flux_theta_global(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	_, _, _, _, rho_theta_ll, _ = u_ll
	_, _, _, _, rho_theta_rr, _ = u_rr

	gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr, equations)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v3_avg = 0.5f0 * (v3_ll + v3_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	theta_mean = ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	if orientation == 1
		f5 = gammamean * v1_avg
		f1 = theta_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg
	elseif orientation == 2
		f5 = gammamean * v2_avg
		f1 = theta_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 * v3_avg
	else
		f5 = gammamean * v3_avg
		f1 = theta_mean * v3_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg
		f4 = f1 * v3_avg + p_avg
	end

	return SVector(f1, f2, f3, f4, f5, zero(eltype(u)))
end

@inline function prim2cons(prim,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, v1, v2, v3, p, phi = prim
	rho_v1 = rho * v1
	rho_v2 = rho * v2
	rho_v3 = rho * v3
	rho_theta = (p / equations.K)^(1 / equations.gamma)
	return SVector(rho, rho_v1, rho_v2, rho_v3, rho_theta, phi)
end

@inline function cons2prim(u,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, rho_v1, rho_v2, rho_v3, rho_theta, phi = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v3 = rho_v3 / rho
	# if rho_theta < 0 
	#     @show rho_theta
	# end
	p = equations.K * rho_theta^equations.gamma
	#p = equations.K * exp(equations.gamma * log(rho_theta))
	# @show "Disp time"
	# elapsed_time_1 = @elapsed p = equations.K*rho_theta^equations.gamma
	# elapsed_time_2 = @elapsed p = equations.p_0*exp(equations.R*rho_theta/equations.p_0 * log(equations.gamma))
	# @show elapsed_time_1
	# @show elapsed_time_2

	# @btime p = equations.K*rho_theta^equations.gamma
	# @btime p = equations.p_0*exp(equations.R*rho_theta/equations.p_0 * log(equations.gamma))
	return SVector(rho, v1, v2, v3, p, phi)
end

@inline function cons2cons(u,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	return u
end

@inline function cons2entropy(u,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, rho_v1, rho_v2, rho_v3, rho_theta, phi = u

	w1 = -0.5f0 * (rho_v1^2 + rho_v2^2 + rho_v3^2) / rho^2
	w2 = rho_v1 / rho
	w3 = rho_v2 / rho
	w4 = rho_v3 / rho
	w5 = equations.gamma * equations.inv_gamma_minus_one * equations.K * (rho_theta)^(equations.gamma - 1)

	return SVector(w1, w2, w3, w4, w5, phi)
end

@inline function cons2entropy2(u, equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, rho_v1, rho_v2, rho_v3, rho_theta, phi = u

	w1 = log(equations.K * (rho_theta / rho)^equations.gamma) - equations.gamma
	w2 = 0.0
	w3 = 0.0
	w4 = 0.0
	w5 = rho / rho_theta * equations.gamma

	return SVector(w1, w2, w3, w4, w5, phi)
end

@inline function entropy_math(cons,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	# Mathematical entropy
	p = equations.K * cons[5]^equations.gamma
	U = (p / (equations.gamma - 1) + 0.5f0 * (cons[2]^2 + cons[3]^2 + cons[4]^2) / (cons[1]))

	return U
end


@inline function entropy_phys(cons, equations::CompressibleEulerPotentialTemperatureEquations3DNC)

	p = equations.K * cons[5]^equations.gamma
	# Thermodynamic entropy
	s = log(p) - equations.gamma * log(cons[1])
	S = -s * cons[1] / (equations.gamma - 1.0)
	return S
end

# Default entropy is the mathematical entropy
@inline function entropy(cons,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	entropy_math(cons, equations)
end

@inline function energy_total(cons,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	entropy(cons, equations)
end

@inline function energy_kinetic(cons,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	return 0.5f0 * (cons[2]^2 + cons[3]^2 + cons[4]^2) / (cons[1])
end

@inline function max_abs_speeds(u, equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, v1, v2, v3, p, _ = cons2prim(u, equations)
	c = sqrt(equations.gamma * p / rho)

	return abs(v1) + c, abs(v2) + c, abs(v3) + c
end

@inline function density_pressure(u, equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho, rho_v1, rho_v2, rho_v3, rho_theta, _ = u
	rho_times_p = rho * equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	return rho_times_p
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	# Get the velocity value in the appropriate direction
	if orientation == 1
		v_ll = v1_ll
		v_rr = v1_rr
	elseif orientation == 2
		v_ll = v2_ll
		v_rr = v2_rr
	else # orientation == 3
		v_ll = v3_ll
		v_rr = v3_rr
	end
	# Calculate sound speeds
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	λ_max = max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr)
end

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	rho_ll, v1_ll, v2_ll, v3_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, v3_rr, p_rr, _ = cons2prim(u_rr, equations)

	# Calculate normal velocities and sound speed
	# left
	v_ll = (v1_ll * normal_direction[1]
			+ v2_ll * normal_direction[2]
			+ v3_ll * normal_direction[3])
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	# right
	v_rr = (v1_rr * normal_direction[1]
			+ v2_rr * normal_direction[2]
			+ v3_rr * normal_direction[3])
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * Trixi.norm(normal_direction)
end

@inline function pressure(cons, equations::CompressibleEulerPotentialTemperatureEquations3DNC)
	_, _, _, _, p, _ = cons2prim(cons, equations)
	return p

end

end # @muladd