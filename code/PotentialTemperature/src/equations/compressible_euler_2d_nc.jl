# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations, norm
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, flux_ranocha, have_nonconservative_terms, boundary_condition_slip_wall, flux_shima_etal, flux_kennedy_gruber, flux_chandrashekar, flux_hllc
@muladd begin
#! format: noindent

@doc raw"""
	CompressibleEulerEquations2D(gamma)

The compressible Euler equations
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho v_2 \\ \rho e
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
 \rho v_1 \\ \rho v_1^2 + p \\ \rho v_1 v_2 \\ (\rho e +p) v_1
\end{pmatrix}
+
\frac{\partial}{\partial y}
\begin{pmatrix}
\rho v_2 \\ \rho v_1 v_2 \\ \rho v_2^2 + p \\ (\rho e +p) v_2
\end{pmatrix}
=
\begin{pmatrix}
0 \\ 0 \\ 0 \\ 0
\end{pmatrix}
```
for an ideal gas with ratio of specific heats `gamma`
in two space dimensions.
Here, ``\rho`` is the density, ``v_1``, ``v_2`` the velocities, ``e`` the specific total energy **rather than** specific internal energy, and
```math
p = (\gamma - 1) \left( \rho e - \frac{1}{2} \rho (v_1^2+v_2^2) \right)
```
the pressure.
"""
struct CompressibleEulerEquations2DNC{RealT <: Real} <:
	   AbstractCompressibleEulerEquations{2, 5}
	gamma::RealT               # ratio of specific heats
	inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications

	function CompressibleEulerEquations2DNC(gamma)
		γ, inv_gamma_minus_one = promote(gamma, inv(gamma - 1))
		new{typeof(γ)}(γ, inv_gamma_minus_one)
	end
end

function varnames(::typeof(cons2cons), ::CompressibleEulerEquations2DNC)
	("rho", "rho_v1", "rho_v2", "rho_e", "phi")
end
varnames(::typeof(cons2prim), ::CompressibleEulerEquations2DNC) = ("rho", "v1", "v2", "p", "phi")

have_nonconservative_terms(::CompressibleEulerEquations2DNC) = Trixi.True()

@inline function flux_nonconservative_gravity_log(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	# Pull the necessary left and right state information
	rho_ll, rho_v1_ll, rho_v2_ll, _, _ = u_ll
	rho_rr, rho_v1_rr, rho_v2_rr, _, _ = u_rr
	rho_avg = ln_mean(rho_ll, rho_rr)
	v1_ll = rho_v1_ll / rho_ll
	v1_rr = rho_v1_rr / rho_rr
	v2_ll = rho_v2_ll / rho_ll
	v2_rr = rho_v2_rr / rho_rr
	v2_avg = 0.5f0 * (v2_rr + v2_ll)
	v1_avg = 0.5f0 * (v1_rr + v1_ll)
	phi_jump = u_rr[5] - u_ll[5]
	RealT = eltype(rho_ll)
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	return SVector(zero(RealT),
		normal_direction[1] * RealT(9.81) * rho_avg * phi_jump,
		normal_direction[2] * RealT(9.81) * rho_avg * phi_jump,
		normal_direction[1] * RealT(9.81) * rho_avg * phi_jump * v1_avg + normal_direction[2] * RealT(9.81) * rho_avg * phi_jump * v2_avg, zero(RealT))
end

@inline function flux_nonconservative_gravity_am(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	RealT = eltype(u_ll)
	# Pull the necessary left and right state information
	rho_ll, rho_v1_ll, rho_v2_ll, _, _ = u_ll
	rho_rr, rho_v1_rr, rho_v2_rr, _, _ = u_rr
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v2_ll = rho_v2_ll / rho_ll
	v2_rr = rho_v2_rr / rho_rr
	v2_avg = 0.5f0 * (v2_rr + v2_ll)
	v1_avg = 0.5f0 * (rho_v1_rr / rho_rr + rho_v1_ll / rho_ll)
	phi_jump = u_rr[5] - u_ll[5]
	return SVector(0,
		normal_direction[1] * RealT(9.81) * rho_avg * phi_jump,
		normal_direction[2] * RealT(9.81) * rho_avg * phi_jump,
		normal_direction[1] * RealT(9.81) * rho_avg * phi_jump * v1_avg + normal_direction[2] * RealT(9.81) * rho_avg * phi_jump * v2_avg, 0)
end

@inline function flux_nonconservative_gravity_gamma(u_ll, u_rr,
	normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	# Pull the necessary left and right state information
	RealT = eltype(u_ll)
	rho_ll, rho_v1_ll, rho_v2_ll, _, _ = u_ll
	rho_rr, rho_v1_rr, rho_v2_rr, _, _ = u_rr
	rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
	v2_ll = rho_v2_ll / rho_ll
	v2_rr = rho_v2_rr / rho_rr
	v2_avg = 0.5f0 * (v2_rr + v2_ll)
	v1_avg = 0.5f0 * (rho_v1_rr / rho_rr + rho_v1_ll / rho_ll)
	phi_jump = u_rr[5] - u_ll[5]
	# Bottom gradient nonconservative term: (0, g h b_x, g h b_y, 0)
	return SVector(0,
		normal_direction[1] * 9.81 * rho_avg * phi_jump,
		normal_direction[2] * 9.81 * rho_avg * phi_jump,
		normal_direction[1] * RealT(9.81) * rho_avg * phi_jump * v1_avg + normal_direction[2] * RealT(9.81) * rho_avg * phi_jump * v2_avg, 0)
end

"""
	boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
								 surface_flux_function, equations::CompressibleEulerEquations2DNC)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, orientation,
	direction, x, t,
	surface_flux_function,
	equations::CompressibleEulerEquations2DNC)
	# get the appropriate normal vector from the orientation
	RealT = eltype(u_inner)
	if orientation == 1
		normal_direction = SVector(one(RealT), zero(RealT))
	else # orientation == 2
		normal_direction = SVector(zero(RealT), one(RealT))
	end

	# compute and return the flux using `boundary_condition_slip_wall` routine above
	return boundary_condition_slip_wall(u_inner, normal_direction, direction,
		x, t, surface_flux_function, equations)
end

"""
	boundary_condition_slip_wall(u_inner, normal_direction, direction, x, t,
								 surface_flux_function, equations::CompressibleEulerEquations2DNC)

Should be used together with [`StructuredMesh`](@ref).
"""
@inline function boundary_condition_slip_wall(u_inner, normal_direction::AbstractVector,
	direction, x, t,
	surface_flux_function,
	equations::CompressibleEulerEquations2DNC)
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


@inline function boundary_condition_slip_wall_2(u_inner, normal_direction::AbstractVector,
	x, t,
	surface_flux_functions,
	equations::CompressibleEulerEquations2DNC)
	surface_flux_function, nonconservative_flux_function = surface_flux_functions
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
	noncons_flux = nonconservative_flux_function(u_inner, u_boundary, normal_direction,
		equations)

	return flux, noncons_flux
end

"""
boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
surface_flux_function, equations::ShallowWaterEquations2D)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall_2(u_inner, orientation,
	direction, x, t,
	surface_flux_functions,
	equations::CompressibleEulerEquations2DNC)
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
		noncons_flux = nonconservative_flux_function(u_inner, u_boundary, orientation,
			equations)
	else # u_boundary is "left" of boundary, u_inner is "right" of boundary
		flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
		noncons_flux = nonconservative_flux_function(u_boundary, u_inner, orientation,
			equations)
	end

	return flux, noncons_flux
end


# Calculate 2D flux for a single point
@inline function flux(u, orientation::Integer, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e, _ = u
	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2))
	if orientation == 1
		f1 = rho_v1
		f2 = rho_v1 * v1 + p
		f3 = rho_v1 * v2
		f4 = (rho_e + p) * v1
	else
		f1 = rho_v2
		f2 = rho_v2 * v1
		f3 = rho_v2 * v2 + p
		f4 = (rho_e + p) * v2
	end
	return SVector(f1, f2, f3, f4, 0)
end

# Calculate 2D flux for a single point in the normal direction
# Note, this directional vector is not normalized
@inline function flux(u, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	rho_e = u[4]
	rho, v1, v2, p, _ = cons2prim(u, equations)

	v_normal = v1 * normal_direction[1] + v2 * normal_direction[2]
	rho_v_normal = rho * v_normal
	f1 = rho_v_normal
	f2 = rho_v_normal * v1 + p * normal_direction[1]
	f3 = rho_v_normal * v2 + p * normal_direction[2]
	f4 = (rho_e + p) * v_normal
	return SVector(f1, f2, f3, f4, 0)
end

"""
	flux_shima_etal(u_ll, u_rr, orientation_or_normal_direction,
					equations::CompressibleEulerEquations2DNC)

This flux is is a modification of the original kinetic energy preserving two-point flux by
- Yuichi Kuya, Kosuke Totani and Soshi Kawai (2018)
  Kinetic energy and entropy preserving schemes for compressible flows
  by split convective forms
  [DOI: 10.1016/j.jcp.2018.08.058](https://doi.org/10.1016/j.jcp.2018.08.058)

The modification is in the energy flux to guarantee pressure equilibrium and was developed by
- Nao Shima, Yuichi Kuya, Yoshiharu Tamaki, Soshi Kawai (JCP 2020)
  Preventing spurious pressure oscillations in split convective form discretizations for
  compressible flows
  [DOI: 10.1016/j.jcp.2020.110060](https://doi.org/10.1016/j.jcp.2020.110060)
"""
@inline function flux_shima_etal(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	kin_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

	# Calculate fluxes depending on orientation
	if orientation == 1
		pv1_avg = 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)
		f1 = rho_avg * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg
	else
		pv2_avg = 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)
		f1 = rho_avg * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = p_avg * v2_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv2_avg
	end

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_shima_etal(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v_dot_n_avg = 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_avg * v_dot_n_avg
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = (f1 * velocity_square_avg +
		  p_avg * v_dot_n_avg * equations.inv_gamma_minus_one
		  + 0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll))

	return SVector(f1, f2, f3, f4, 0)
end

"""
	flux_kennedy_gruber(u_ll, u_rr, orientation_or_normal_direction,
						equations::CompressibleEulerEquations2DNC)

Kinetic energy preserving two-point flux by
- Kennedy and Gruber (2008)
  Reduced aliasing formulations of the convective terms within the
  Navier-Stokes equations for a compressible fluid
  [DOI: 10.1016/j.jcp.2007.09.020](https://doi.org/10.1016/j.jcp.2007.09.020)
"""
@inline function flux_kennedy_gruber(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state
	rho_e_ll = u_ll[4]
	rho_e_rr = u_rr[4]
	rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

	# Calculate fluxes depending on orientation
	if orientation == 1
		f1 = rho_avg * v1_avg
		f2 = rho_avg * v1_avg * v1_avg + p_avg
		f3 = rho_avg * v1_avg * v2_avg
		f4 = (rho_avg * e_avg + p_avg) * v1_avg
	else
		f1 = rho_avg * v2_avg
		f2 = rho_avg * v2_avg * v1_avg
		f3 = rho_avg * v2_avg * v2_avg + p_avg
		f4 = (rho_avg * e_avg + p_avg) * v2_avg
	end

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_kennedy_gruber(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state
	rho_e_ll = u_ll[4]
	rho_e_rr = u_rr[4]
	rho_ll, v1_ll, v2_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr = cons2prim(u_rr, equations)

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	v_dot_n_avg = v1_avg * normal_direction[1] + v2_avg * normal_direction[2]
	p_avg = 0.5f0 * (p_ll + p_rr)
	e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_avg * v_dot_n_avg
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = f1 * e_avg + p_avg * v_dot_n_avg

	return SVector(f1, f2, f3, f4, 0)
end

"""
	flux_chandrashekar(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquations2DNC)

Entropy conserving two-point flux by
- Chandrashekar (2013)
  Kinetic Energy Preserving and Entropy Stable Finite Volume Schemes
  for Compressible Euler and Navier-Stokes Equations
  [DOI: 10.4208/cicp.170712.010313a](https://doi.org/10.4208/cicp.170712.010313a)
"""
@inline function flux_chandrashekar(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	beta_ll = 0.5f0 * rho_ll / p_ll
	beta_rr = 0.5f0 * rho_rr / p_rr
	specific_kin_ll = 0.5f0 * (v1_ll^2 + v2_ll^2)
	specific_kin_rr = 0.5f0 * (v1_rr^2 + v2_rr^2)

	# Compute the necessary mean values
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	rho_mean = ln_mean(rho_ll, rho_rr)
	beta_mean = ln_mean(beta_ll, beta_rr)
	beta_avg = 0.5f0 * (beta_ll + beta_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_mean = 0.5f0 * rho_avg / beta_avg
	velocity_square_avg = specific_kin_ll + specific_kin_rr

	# Calculate fluxes depending on orientation
	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_mean
		f3 = f1 * v2_avg
		f4 = f1 * 0.5f0 *
			 (1 / (equations.gamma - 1) / beta_mean - velocity_square_avg) +
			 f2 * v1_avg + f3 * v2_avg
	else
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_mean
		f4 = f1 * 0.5f0 *
			 (1 / (equations.gamma - 1) / beta_mean - velocity_square_avg) +
			 f2 * v1_avg + f3 * v2_avg
	end

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_chandrashekar(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]
	beta_ll = 0.5f0 * rho_ll / p_ll
	beta_rr = 0.5f0 * rho_rr / p_rr
	specific_kin_ll = 0.5f0 * (v1_ll^2 + v2_ll^2)
	specific_kin_rr = 0.5f0 * (v1_rr^2 + v2_rr^2)

	# Compute the necessary mean values
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	rho_mean = ln_mean(rho_ll, rho_rr)
	beta_mean = ln_mean(beta_ll, beta_rr)
	beta_avg = 0.5f0 * (beta_ll + beta_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_mean = 0.5f0 * rho_avg / beta_avg
	velocity_square_avg = specific_kin_ll + specific_kin_rr

	# Multiply with average of normal velocities
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_mean * normal_direction[1]
	f3 = f1 * v2_avg + p_mean * normal_direction[2]
	f4 = f1 * 0.5f0 * (1 / (equations.gamma - 1) / beta_mean - velocity_square_avg) +
		 f2 * v1_avg + f3 * v2_avg

	return SVector(f1, f2, f3, f4, 0)
end

"""
	flux_ranocha(u_ll, u_rr, orientation_or_normal_direction,
				 equations::CompressibleEulerEquations2DNC)

Entropy conserving and kinetic energy preserving two-point flux by
- Hendrik Ranocha (2018)
  Generalised Summation-by-Parts Operators and Entropy Stability of Numerical Methods
  for Hyperbolic Balance Laws
  [PhD thesis, TU Braunschweig](https://cuvillier.de/en/shop/publications/7743)
See also
- Hendrik Ranocha (2020)
  Entropy Conserving and Kinetic Energy Preserving Numerical Methods for
  the Euler Equations Using Summation-by-Parts Operators
  [Proceedings of ICOSAHOM 2018](https://doi.org/10.1007/978-3-030-39647-3_42)
"""
@inline function flux_ranocha(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)

	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	# Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
	# in exact arithmetic since
	#     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
	#   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
	inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

	# Calculate fluxes depending on orientation
	if orientation == 1
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = f1 * v2_avg
		f4 = f1 *
			 (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
			 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)
	else
		f1 = rho_mean * v2_avg
		f2 = f1 * v1_avg
		f3 = f1 * v2_avg + p_avg
		f4 = f1 *
			 (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
			 0.5f0 * (p_ll * v2_rr + p_rr * v2_ll)
	end

	return SVector(f1, f2, f3, f4, 0)
end

@inline function flux_ranocha(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	# Unpack left and right state

	rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)
	v_dot_n_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_dot_n_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	# Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
	# in exact arithmetic since
	#     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
	#   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
	inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	v2_avg = 0.5f0 * (v2_ll + v2_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	velocity_square_avg = 0.5f0 * (v1_ll * v1_rr + v2_ll * v2_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * 0.5f0 * (v_dot_n_ll + v_dot_n_rr)
	f2 = f1 * v1_avg + p_avg * normal_direction[1]
	f3 = f1 * v2_avg + p_avg * normal_direction[2]
	f4 = (f1 * (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one)
		  +
		  0.5f0 * (p_ll * v_dot_n_rr + p_rr * v_dot_n_ll))

	return SVector(f1, f2, f3, f4, 0)
end


@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
	c = flux_lmars.speed_of_sound

	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)

	if orientation == 1
		v_ll = v1_ll
		v_rr = v1_rr
	else # orientation == 2
		v_ll = v2_ll
		v_rr = v2_rr
	end

	rho = 0.5f0 * (rho_ll + rho_rr)
	p = 0.5f0 * (p_ll + p_rr) - 0.5f0 * c * rho * (v_rr - v_ll)
	v = 0.5f0 * (v_ll + v_rr) - 0.5f0 / (c * rho) * (p_rr - p_ll)

	# We treat the energy term analogous to the potential temperature term in the paper by
	# Chen et al., i.e. we use p_ll and p_rr, and not p
	if v >= 0
		f1, f2, f3, f4 = v * u_ll
		f4 = f4 + p_ll * v
	else
		f1, f2, f3, f4 = v * u_rr
		f4 = f4 + p_rr * v
	end

	if orientation == 1
		f2 = f2 + p
	else # orientation == 2
		f3 = f3 + p
	end

	return SVector(f1, f2, f3, f4, 0)
end

@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	c = flux_lmars.speed_of_sound

	# Unpack left and right state
	rho_ll, v1_ll, v2_ll, p_ll, phi_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, phi_rr = cons2prim(u_rr, equations)

	v_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	# Note that this is the same as computing v_ll and v_rr with a normalized normal vector
	# and then multiplying v by `norm_` again, but this version is slightly faster.
	norm_ = norm(normal_direction)

	rho = 0.5f0 * (rho_ll + rho_rr)
	p = 0.5f0 * (p_ll + p_rr) - 0.5f0 * c * rho * (v_rr - v_ll) / norm_
	v = 0.5f0 * (v_ll + v_rr) - 0.5f0 / (c * rho) * (p_rr - p_ll) * norm_

	# We treat the energy term analogous to the potential temperature term in the paper by
	# Chen et al., i.e. we use p_ll and p_rr, and not p
	if v >= 0
		f1, f2, f3, f4 = u_ll * v
		f4 = f4 + p_ll * v
	else
		f1, f2, f3, f4 = u_rr * v
		f4 = f4 + p_rr * v
	end

	return SVector(f1,
		f2 + p * normal_direction[1],
		f3 + p * normal_direction[2],
		f4, 0.0)
end

# Calculate maximum wave speed for local Lax-Friedrichs-type dissipation as the
# maximum velocity magnitude plus the maximum speed of sound
@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
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

@inline function max_abs_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	# Calculate normal velocities and sound speed
	# left
	v_ll = (v1_ll * normal_direction[1]
			+
			v2_ll * normal_direction[2])
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	# right
	v_rr = (v1_rr * normal_direction[1]
			+
			v2_rr * normal_direction[2])
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	return max(abs(v_ll), abs(v_rr)) + max(c_ll, c_rr) * norm(normal_direction)
end

# Calculate estimate for minimum and maximum wave speeds for HLL-type fluxes
@inline function min_max_speed_naive(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	if orientation == 1 # x-direction
		λ_min = v1_ll - sqrt(equations.gamma * p_ll / rho_ll)
		λ_max = v1_rr + sqrt(equations.gamma * p_rr / rho_rr)
	else # y-direction
		λ_min = v2_ll - sqrt(equations.gamma * p_ll / rho_ll)
		λ_max = v2_rr + sqrt(equations.gamma * p_rr / rho_rr)
	end

	return λ_min, λ_max
end

@inline function min_max_speed_naive(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	v_normal_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_normal_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	norm_ = norm(normal_direction)
	# The v_normals are already scaled by the norm
	λ_min = v_normal_ll - sqrt(equations.gamma * p_ll / rho_ll) * norm_
	λ_max = v_normal_rr + sqrt(equations.gamma * p_rr / rho_rr) * norm_

	return λ_min, λ_max
end

# More refined estimates for minimum and maximum wave speeds for HLL-type fluxes
@inline function min_max_speed_davis(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations2DNC)
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	if orientation == 1 # x-direction
		λ_min = min(v1_ll - c_ll, v1_rr - c_rr)
		λ_max = max(v1_ll + c_ll, v1_rr + c_rr)
	else # y-direction
		λ_min = min(v2_ll - c_ll, v2_rr - c_rr)
		λ_max = max(v2_ll + c_ll, v2_rr + c_rr)
	end

	return λ_min, λ_max
end

# More refined estimates for minimum and maximum wave speeds for HLL-type fluxes
@inline function min_max_speed_davis(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations2DNC)
	rho_ll, v1_ll, v2_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, v2_rr, p_rr, _ = cons2prim(u_rr, equations)

	norm_ = norm(normal_direction)

	c_ll = sqrt(equations.gamma * p_ll / rho_ll) * norm_
	c_rr = sqrt(equations.gamma * p_rr / rho_rr) * norm_

	v_normal_ll = v1_ll * normal_direction[1] + v2_ll * normal_direction[2]
	v_normal_rr = v1_rr * normal_direction[1] + v2_rr * normal_direction[2]

	# The v_normals are already scaled by the norm
	λ_min = min(v_normal_ll - c_ll, v_normal_rr - c_rr)
	λ_max = max(v_normal_ll + c_ll, v_normal_rr + c_rr)

	return λ_min, λ_max
end

# Called inside `FluxRotated` in `numerical_fluxes.jl` so the direction
# has been normalized prior to this rotation of the state vector
@inline function rotate_to_x(u, normal_vector, equations::CompressibleEulerEquations2DNC)
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

# Called inside `FluxRotated` in `numerical_fluxes.jl` so the direction
# has been normalized prior to this back-rotation of the state vector
@inline function rotate_from_x(u, normal_vector,
	equations::CompressibleEulerEquations2DNC)
	# cos and sin of the angle between the x-axis and the normalized normal_vector are
	# the normalized vector's x and y coordinates respectively (see unit circle).
	c = normal_vector[1]
	s = normal_vector[2]

	# Apply the 2D back-rotation matrix with normal and tangent directions of the form
	# [ 1    0    0   0;
	#   0   n_1  t_1  0;
	#   0   n_2  t_2  0;
	#   0    0    0   1 ]
	# where t_1 = -n_2 and t_2 = n_1

	return SVector(u[1],
		c * u[2] - s * u[3],
		s * u[2] + c * u[3],
		u[4], u[5])
end

@inline function max_abs_speeds(u, equations::CompressibleEulerEquations2DNC)
	rho, v1, v2, p, phi = cons2prim(u, equations)
	c = sqrt(equations.gamma * p / rho)

	return abs(v1) + c, abs(v2) + c
end

# Convert conservative variables to primitive
@inline function cons2prim(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e, phi = u

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2))

	return SVector(rho, v1, v2, p, phi)
end

# Convert conservative variables to entropy
@inline function cons2entropy(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e, _ = u

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v_square = v1^2 + v2^2
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * v_square)
	s = log(p) - equations.gamma * log(rho)
	rho_p = rho / p

	w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
		 0.5f0 * rho_p * v_square
	w2 = rho_p * v1
	w3 = rho_p * v2
	w4 = -rho_p

	return SVector(w1, w2, w3, w4, 0)
end

# Transformation from conservative variables u to entropy vector ds_0/du,
# using the modified specific entropy of Guermond et al. (2019): s_0 = p * rho^(-gamma) / (gamma-1).
# Note: This is *not* the "conventional" specific entropy s = ln(p / rho^(gamma)).
@inline function cons2entropy_guermond_etal(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e, phi = u

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v_square = v1^2 + v2^2
	inv_rho_gammap1 = (1 / rho)^(equations.gamma + 1)

	# The derivative vector for the modified specific entropy of Guermond et al.
	w1 = inv_rho_gammap1 *
		 (0.5f0 * rho * (equations.gamma + 1) * v_square - equations.gamma * rho_e)
	w2 = -rho_v1 * inv_rho_gammap1
	w3 = -rho_v2 * inv_rho_gammap1
	w4 = (1 / rho)^equations.gamma

	return SVector(w1, w2, w3, w4, 0)
end

@inline function entropy2cons(w, equations::CompressibleEulerEquations2DNC)
	# See Hughes, Franca, Mallet (1986) A new finite element formulation for CFD
	# [DOI: 10.1016/0045-7825(86)90127-1](https://doi.org/10.1016/0045-7825(86)90127-1)
	@unpack gamma = equations

	# convert to entropy `-rho * s` used by Hughes, France, Mallet (1986)
	# instead of `-rho * s / (gamma - 1)`
	V1, V2, V3, V5, phi = w .* (gamma - 1)

	# s = specific entropy, eq. (53)
	s = gamma - V1 + (V2^2 + V3^2) / (2 * V5)

	# eq. (52)
	rho_iota = ((gamma - 1) / (-V5)^gamma)^(equations.inv_gamma_minus_one) *
			   exp(-s * equations.inv_gamma_minus_one)

	# eq. (51)
	rho = -rho_iota * V5
	rho_v1 = rho_iota * V2
	rho_v2 = rho_iota * V3
	rho_e = rho_iota * (1 - (V2^2 + V3^2) / (2 * V5))
	return SVector(rho, rho_v1, rho_v2, rho_e, phi)
end

# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::CompressibleEulerEquations2DNC)
	rho, v1, v2, p, phi = prim
	rho_v1 = rho * v1
	rho_v2 = rho * v2
	rho_e = p * equations.inv_gamma_minus_one + 0.5f0 * (rho_v1 * v1 + rho_v2 * v2)
	return SVector(rho, rho_v1, rho_v2, rho_e, phi)
end

@inline function density(u, equations::CompressibleEulerEquations2DNC)
	rho = u[1]
	return rho
end

@inline function velocity(u, equations::CompressibleEulerEquations2DNC)
	rho = u[1]
	v1 = u[2] / rho
	v2 = u[3] / rho
	return SVector(v1, v2)
end

@inline function velocity(u, orientation::Int, equations::CompressibleEulerEquations2DNC)
	rho = u[1]
	v = u[orientation+1] / rho
	return v
end

@inline function pressure(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e, _ = u
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1^2 + rho_v2^2) / rho)
	return p
end

# Transformation from conservative variables u to d(p)/d(u)
@inline function gradient_conservative(::typeof(pressure),
	u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e, _ = u

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v_square = v1^2 + v2^2

	return (equations.gamma - 1) * SVector(0.5f0 * v_square, -v1, -v2, 1)
end

@inline function density_pressure(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e, _ = u
	rho_times_p = (equations.gamma - 1) * (rho * rho_e - 0.5f0 * (rho_v1^2 + rho_v2^2))
	return rho_times_p
end

# Calculates the entropy flux in direction "orientation" and the entropy variables for a state cons
# NOTE: This method seems to work currently (b82534e) but is never used anywhere. Thus it is
# commented here until someone uses it or writes a test for it.
# @inline function cons2entropyvars_and_flux(gamma::Float64, cons, orientation::Int)
#   entropy = MVector{4, Float64}(undef)
#   v = (cons[2] / cons[1] , cons[3] / cons[1])
#   v_square= v[1]*v[1]+v[2]*v[2]
#   p = (gamma - 1) * (cons[4] - 1/2 * (cons[2] * v[1] + cons[3] * v[2]))
#   rho_p = cons[1] / p
#   # thermodynamic entropy
#   s = log(p) - gamma*log(cons[1])
#   # mathematical entropy
#   S = - s*cons[1]/(gamma-1)
#   # entropy variables
#   entropy[1] = (gamma - s)/(gamma-1) - 0.5*rho_p*v_square
#   entropy[2] = rho_p*v[1]
#   entropy[3] = rho_p*v[2]
#   entropy[4] = -rho_p
#   # entropy flux
#   entropy_flux = S*v[orientation]
#   return entropy, entropy_flux
# end

# Calculate thermodynamic entropy for a conservative state `cons`
@inline function entropy_thermodynamic(cons, equations::CompressibleEulerEquations2DNC)
	# Pressure
	p = (equations.gamma - 1) * (cons[4] - 0.5f0 * (cons[2]^2 + cons[3]^2) / cons[1])

	# Thermodynamic entropy
	s = log(p) - equations.gamma * log(cons[1])

	return s
end

# Calculate mathematical entropy for a conservative state `cons`
@inline function entropy_math(cons, equations::CompressibleEulerEquations2DNC)
	# Mathematical entropy
	S = -entropy_thermodynamic(cons, equations) * cons[1] *
		equations.inv_gamma_minus_one

	return S
end

# Transformation from conservative variables u to d(s)/d(u)
@inline function gradient_conservative(::typeof(entropy_math),
	u, equations::CompressibleEulerEquations2DNC)
	return cons2entropy(u, equations)
end


# Default entropy is the mathematical entropy
@inline function entropy(cons, equations::CompressibleEulerEquations2DNC)
	entropy_math(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline energy_total(cons, ::CompressibleEulerEquations2DNC) = cons[4]

# Calculate kinetic energy for a conservative state `cons`
@inline function energy_kinetic(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_e = u
	return (rho_v1^2 + rho_v2^2) / (2 * rho)
end

# Calculate internal energy for a conservative state `cons`
@inline function energy_internal(cons, equations::CompressibleEulerEquations2DNC)
	return energy_total(cons, equations) - energy_kinetic(cons, equations)
end

@inline function well_balanced_v1(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, _ = u
	v1 = rho_v1 / rho
	return abs(v1)
end

@inline function well_balanced_v2(u, equations::CompressibleEulerEquations2DNC)
	rho, rho_v1, rho_v2, rho_theta, _ = u
	v2 = rho_v2 / rho
	return abs(v2)
end

end # @muladd
