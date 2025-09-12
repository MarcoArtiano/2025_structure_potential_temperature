using Trixi
using Trixi: ln_mean, stolarsky_mean
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, boundary_condition_slip_wall

@muladd begin
#! format: noindent
struct CompressibleEulerPotentialTemperatureEquations1DNC{RealT <: Real} <:
	   Trixi.AbstractCompressibleEulerEquations{1, 4}
	p_0::RealT
	c_p::RealT
	c_v::RealT
	g::RealT
	R::RealT
	gamma::RealT
	a::RealT
	inv_gamma_minus_one::RealT
	K::RealT
end

function CompressibleEulerPotentialTemperatureEquations1DNC(; g = 9.81, RealT = Float64)
	p_0 = 100_000.0
	c_p = 1004.0
	c_v = 717.0
	R = c_p - c_v
	gamma = c_p / c_v
	a = 340.0
	inv_gamma_minus_one = inv(gamma - 1)
	K = p_0 * (R / p_0)^gamma

	return CompressibleEulerPotentialTemperatureEquations1DNC{RealT}(p_0, c_p, c_v, g, R, gamma, a, inv_gamma_minus_one, K)
end

function varnames(::typeof(cons2cons), ::CompressibleEulerPotentialTemperatureEquations1DNC)
	("rho", "rho_v1", "rho_theta", "phi")
end

varnames(::typeof(cons2prim), ::CompressibleEulerPotentialTemperatureEquations1DNC) = ("rho", "v1", "p1", "phi")

Trixi.have_nonconservative_terms(::CompressibleEulerPotentialTemperatureEquations1DNC) = Trixi.True()

@inline function flux(u, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho, rho_v1, rho_theta, phi = u
	v1 = rho_v1 / rho
	p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	p = equations.K * exp(log(rho_theta^equations.gamma))

	f1 = rho_v1
	f2 = rho_v1 * v1 + p
	f3 = rho_theta * v1

	return SVector(f1, f2, f3, zero(eltype(u)))
end

@inline function flux_nonconservative_gravity_am(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, phi_ll = u_ll
	rho_rr, _, _, phi_rr = u_rr

	rho_avg = 0.5f0 * (rho_ll + rho_rr)

	jphi = phi_rr - phi_ll

	# Bottom gradient nonconservative term: (0, g h b_x, 0)
	f = SVector(0.0, rho_avg * jphi, 0.0, 0.0)

	return f
end

@inline function flux_nonconservative_gravity_log(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, phi_ll = u_ll
	rho_rr, _, _, phi_rr = u_rr

	rho_avg = ln_mean(rho_ll, rho_rr)
	jphi = phi_rr - phi_ll

	# Bottom gradient nonconservative term: (0, g h b_x, 0)
	f = SVector(0.0, rho_avg * jphi, 0.0, 0.0)

	return f
end

@inline function flux_nonconservative_gravity_gamma(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Pull the necessary left and right state information
	rho_ll, _, _, phi_ll = u_ll
	rho_rr, _, _, phi_rr = u_rr

	#rho_avg = 0.5f0*(rho_ll + rho_rr)
	#rho_avg = ln_mean(rho_ll, rho_rr)
	rho_avg = stolarsky_mean(rho_ll, rho_rr, equations.gamma)
	jphi = phi_rr - phi_ll

	# Bottom gradient nonconservative term: (0, g h b_x, 0)
	f = SVector(0.0, rho_avg * jphi, 0.0, 0.0)

	return f
end

@inline function initial_condition_density_wave(x, t, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	v1 = 0.1
	rho = 1 + 0.98 * sinpi(2 * (x[1] - t * v1))
	rho_v1 = rho * v1
	p = 20
	rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
	return SVector(rho, rho_v1, rho_theta, zero(eltype(u)))
end

@inline function initial_condition_weak_blast_wave(x, t,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# From Hennemann & Gassner JCP paper 2020 (Sec. 6.3)
	# Set up polar coordinates
	inicenter = SVector(0.0)
	x_norm = x[1] - inicenter[1]
	r = abs(x_norm)
	# The following code is equivalent to
	# phi = atan(0.0, x_norm)
	# cos_phi = cos(phi)
	# in 1D but faster
	cos_phi = x_norm > 0 ? one(x_norm) : -one(x_norm)

	# Calculate primitive variables
	rho = r > 0.5 ? 1.0 : 1.1691
	v1 = r > 0.5 ? 0.0 : 0.1882 * cos_phi
	p = r > 0.5 ? 1.0 : 1.245

	return prim2cons(SVector(rho, v1, p, zero(eltype(u))), equations)
end

@inline function source_terms_gravity(u, x, t,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho, _, _, _ = u
	return SVector(zero(eltype(u)), -equations.g * rho,
		zero(eltype(u)), zero(eltype(u)))
end

@inline function boundary_condition_slip_wall_2(u_inner, orientation_or_normal, direction,
	x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)

	# create the "external" boundary solution state
	u_boundary = SVector(u_inner[1],
		-u_inner[2],
		u_inner[3],
		u_inner[4])

	# calculate the boundary flux
	if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
		flux = surface_flux_function(u_inner, u_boundary, orientation_or_normal,
			equations)
	else # u_boundary is "left" of boundary, u_inner is "right" of boundary
		flux = surface_flux_function(u_boundary, u_inner, orientation_or_normal,
			equations)
	end

	return flux
end

@inline function boundary_condition_slip_wall(u_inner, orientation,
	direction, x, t,
	surface_flux_function,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# compute the primitive variables
	rho_local, v_normal, p_local = cons2prim(u_inner, equations)

	if isodd(direction) # flip sign of normal to make it outward pointing
		v_normal *= -1
	end

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
		p_star,
		zero(eltype(u_inner)), zero(eltype(u_inner)))
end


# Low Mach number approximate Riemann solver (LMARS) from
# X. Chen, N. Andronova, B. Van Leer, J. E. Penner, J. P. Boyd, C. Jablonowski, S.
# Lin, A Control-Volume Model of the Compressible Euler Equations with a Vertical Lagrangian
# Coordinate Monthly Weather Review Vol. 141.7, pages 2526–2544, 2013,
# https://journals.ametsoc.org/view/journals/mwre/141/7/mwr-d-12-00129.1.xml.
@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	@unpack a = equations
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, _ = cons2prim(u_rr, equations)

	rho = 0.5f0 * (rho_ll + rho_rr)

	p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v1_rr - v1_ll)
	v_interface = 0.5f0 * (v1_ll + v1_rr) - 1 / (2 * a * rho) * (p_rr - p_ll)

	if (v_interface > 0)
		f1, f2, f3, f4 = u_ll * v_interface
	else
		f1, f2, f3, f4 = u_rr * v_interface
	end

	return SVector(f1,
		f2 + p_interface,
		f3, zero(eltype(u_ll)))
end

## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation
@inline function flux_theta(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, _ = cons2prim(u_rr, equations)
	_, _, rho_theta_ll, _ = u_ll
	_, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = gammamean * v1_avg
	return SVector(f1, f2, f3, zero(eltype(u_ll)))
end

@inline function flux_theta_AM(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, _ = cons2prim(u_rr, equations)
	_, _, rho_theta_ll, _ = u_ll
	_, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	#rho_mean = ln_mean(rho_ll, rho_rr)
	rho_mean = 0.5f0 * (rho_ll + rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = gammamean * v1_avg
	return SVector(f1, f2, f3, zero(eltype(u_ll)))
end

@inline function flux_theta_rhos(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, phi = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, phi = cons2prim(u_rr, equations)
	_, _, rho_theta_ll = u_ll
	_, _, rho_theta_rr = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr) * f1
	return SVector(f1, f2, f3, 0)
end

@inline function flux_theta_rhos_AM(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, phi = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, phi = cons2prim(u_rr, equations)
	_, _, rho_theta_ll = u_ll
	_, _, rho_theta_rr = u_rr
	# Compute the necessary mean values
	rho_mean = (rho_ll + rho_rr) * 0.5f0
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr) * f1
	return SVector(f1, f2, f3, 0)
end

@inline function flux_theta_global(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, phi = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, phi = cons2prim(u_rr, equations)
	_, _, rho_theta_ll = u_ll
	_, _, rho_theta_rr = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f3 = gammamean * v1_avg
	f1 = f3 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
	f2 = f1 * v1_avg + p_avg
	return SVector(f1, f2, f3, 0)
end

# Entropy stable, density and pressure positivity preserving flux
@inline function flux_theta_es(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho_ll, v1_ll, p_ll, _ = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, _ = cons2prim(u_rr, equations)
	_, _, rho_theta_ll, _ = u_ll
	_, _, rho_theta_rr, _ = u_rr
	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	rho_mean = 0.5f0 * (rho_ll + rho_rr)
	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)

	# Calculate fluxes depending on normal_direction
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = gammamean * v1_avg

	lambda = max(abs(v1_ll), abs(v1_rr))

	return SVector(f1 - 0.5 * lambda * (rho_rr - rho_ll), f2 - 0.5 * lambda * (rho_rr * v1_rr - rho_ll * v1_ll), f3 - 0.5 * lambda * (rho_theta_rr - rho_theta_ll), zero(eltype(u_ll)))

end

@inline function prim2cons(prim, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho, v1, p, phi = prim
	rho_v1 = rho * v1
	rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
	return SVector(rho, rho_v1, rho_theta, phi)
end

@inline function cons2prim(u, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho, rho_v1, rho_theta, phi = u
	v1 = rho_v1 / rho
	#  p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
	p = equations.K * (rho_theta)^equations.gamma
	p = equations.K * exp(log(rho_theta^equations.gamma))
	return SVector(rho, v1, p, phi)
end

@inline function cons2cons(u, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	return u
end

@inline function cons2entropy(u, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho, rho_v1, rho_theta, phi = u

	k = equations.gamma
	w1 = -0.5f0 * rho_v1^2 / (rho)^2
	w2 = rho_v1 / rho
	w3 = equations.gamma / (equations.gamma - 1) * k * (rho_theta)^(equations.gamma - 1)

	return SVector(w1, w2, w3, zero(eltype(u)))
end

@inline function cons2entropy2(u, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho, rho_v1, rho_theta, phi = u

	k = equations.p_0 * (equations.R / equations.p_0)^equations.gamma
	w1 = log(k * (rho_theta / rho)^equations.gamma) - equations.gamma
	w2 = 0.0
	w3 = rho / rho_theta * equations.gamma

	return SVector(w1, w2, w3, zero(eltype(u)))
end

@inline function entropy_math(cons, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	# Mathematical entropy
	p = equations.p_0 * (equations.R * cons[3] / equations.p_0)^equations.gamma

	U = (p / (equations.gamma - 1) + 1 / 2 * (cons[2]^2) / (cons[1]))

	return U
end

# Default entropy is the mathematical entropy
@inline function entropy(cons, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	entropy_math(cons, equations)
end

@inline function entropy_phys(cons, equations::CompressibleEulerPotentialTemperatureEquations1DNC)

	p = equations.p_0 * (equations.R * cons[3] / equations.p_0)^equations.gamma
	# Thermodynamic entropy
	s = log(p) - equations.gamma * log(cons[1])
	S = -s * cons[1] / (equations.gamma - 1.0)
	return S
end

@inline function energy_total(cons, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	entropy(cons, equations)
end

@inline function energy_kinetic(cons, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	return 0.5f0 * (cons[2]^2) / (cons[1])
end

@inline function max_abs_speeds(u, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho, v1, p, phi = cons2prim(u, equations)
	c = sqrt(equations.gamma * p / rho)

	return (abs(v1) + c,)
end

@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	rho_ll, v1_ll, p_ll, phi = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, phi = cons2prim(u_rr, equations)

	# Calculate primitive variables and speed of sound
	v_mag_ll = abs(v1_ll)
	c_ll = sqrt(equations.gamma * p_ll / rho_ll)
	v_mag_rr = abs(v1_rr)
	c_rr = sqrt(equations.gamma * p_rr / rho_rr)

	λ_max = max(v_mag_ll, v_mag_rr) + max(c_ll, c_rr)
end

@inline function pressurecompute(cons, equations::CompressibleEulerPotentialTemperatureEquations1DNC)
	_, _, p, _ = cons2prim(cons, equations)
	return p

end
end # @muladd
