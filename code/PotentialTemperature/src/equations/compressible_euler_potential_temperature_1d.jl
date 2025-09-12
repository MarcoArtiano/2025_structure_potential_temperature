using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations, inv_ln_mean
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, boundary_condition_slip_wall, energy_total, energy_kinetic, pressure

@muladd begin

	import Trixi: initial_condition_density_wave, flux, max_abs_speed_naive, varnames, cons2prim,
		prim2cons, cons2entropy, density_pressure, max_abs_speeds

#! format: noindent
	struct CompressibleEulerPotentialTemperatureEquations1D{RealT <: Real} <:
		   AbstractCompressibleEulerEquations{1, 3}
		p_0::RealT
		c_p::RealT
		c_v::RealT
		g::RealT
		R::RealT
		gamma::RealT
		a::RealT
		inv_gamma_minus_one::RealT
		K::RealT
		stolarsky_factor::RealT
	end

	function CompressibleEulerPotentialTemperatureEquations1D(; g = 9.81, RealT = Float64)
		p_0 = 100_000.0
		c_p = 1004.0
		c_v = 717.0
		R = c_p - c_v
		gamma = c_p / c_v
		a = 340.0
		inv_gamma_minus_one = inv(gamma - 1)
		K = p_0 * (R / p_0)^gamma
		stolarsky_factor = (gamma - 1.0)/gamma
		return CompressibleEulerPotentialTemperatureEquations1D{RealT}(p_0, c_p, c_v, g, R, gamma, a, inv_gamma_minus_one, K, stolarsky_factor)
	end

	function varnames(::typeof(cons2cons), ::CompressibleEulerPotentialTemperatureEquations1D)
		("rho", "rho_v1", "rho_theta")
	end

	varnames(::typeof(cons2prim), ::CompressibleEulerPotentialTemperatureEquations1D) = ("rho", "v1", "p1")

	@inline function flux(u, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho, rho_v1, rho_theta = u
		v1 = rho_v1 / rho
		p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
		p = equations.K * exp(log(rho_theta^equations.gamma))

		f1 = rho_v1
		f2 = rho_v1 * v1 + p
		f3 = rho_theta * v1

		return SVector(f1, f2, f3)
	end

	@inline function initial_condition_density_wave(x, t, equations::CompressibleEulerPotentialTemperatureEquations1D)
		v1 = 0.1
		rho = 1 + 0.98 * sinpi(2 * (x[1] - t * v1))
		rho_v1 = rho * v1
		p = 20
		rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
		return SVector(rho, rho_v1, rho_theta)
	end

	@inline function initial_condition_weak_blast_wave(x, t,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
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

		return prim2cons(SVector(rho, v1, p), equations)
	end

	@inline function source_terms_gravity(u, x, t,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho, _, _ = u
		return SVector(zero(eltype(u)), -equations.g * rho,
			zero(eltype(u)))
	end

	@inline function boundary_condition_slip_wall(u_inner, orientation,
		direction, x, t,
		surface_flux_function,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
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
			zero(eltype(u_inner)))
	end

	@inline function (flux_lmars::FluxLMARS)(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		a = flux_lmars.speed_of_sound
		rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
		rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

		rho = 0.5f0 * (rho_ll + rho_rr)

		p_interface = 0.5f0 * (p_ll + p_rr) - 0.5f0 * a * rho * (v1_rr - v1_ll)
		v_interface = 0.5f0 * (v1_ll + v1_rr) - 1 / (2 * a * rho) * (p_rr - p_ll)

		if (v_interface > 0)
			f1, f2, f3 = u_ll * v_interface
		else
			f1, f2, f3 = u_rr * v_interface
		end

		return SVector(f1,
			f2 + p_interface,
			f3)
	end

	## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation
	@inline function flux_theta(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		# Unpack left and right state
		rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
		rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
		_, _, rho_theta_ll = u_ll
		_, _, rho_theta_rr = u_rr
		# Compute the necessary mean values
		rho_mean = ln_mean(rho_ll, rho_rr)
		#gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)
		gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)

		v1_avg = 0.5f0 * (v1_ll + v1_rr)
		p_avg = 0.5f0 * (p_ll + p_rr)

		# Calculate fluxes depending on normal_direction
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = gammamean * v1_avg
		return SVector(f1, f2, f3)
	end

	# Entropy stable, density and pressure positivity preserving flux
	@inline function flux_theta_es(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		_, v1_ll, _ = cons2prim(u_ll, equations)
		_, v1_rr, _ = cons2prim(u_rr, equations)

		f_ec = flux_theta(u_ll, u_rr, orientation, equations)

		lambda = max(abs(v1_ll), abs(v1_rr))

		return f_ec - 0.5f0 * lambda * (u_rr - u_ll)
	end

	@inline function flux_theta_AM(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		# Unpack left and right state
		rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
		rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
		_, _, rho_theta_ll = u_ll
		_, _, rho_theta_rr = u_rr
		# Compute the necessary mean values
		rho_mean = 0.5f0 * (rho_ll + rho_rr)
		gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)

		#	gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

		v1_avg = 0.5f0 * (v1_ll + v1_rr)
		p_avg = 0.5f0 * (p_ll + p_rr)

		# Calculate fluxes depending on normal_direction
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = gammamean * v1_avg
		return SVector(f1, f2, f3)

	end

	@inline function flux_theta_AM_es(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		_, v1_ll, _ = cons2prim(u_ll, equations)
		_, v1_rr, _ = cons2prim(u_rr, equations)

		f_ec = flux_theta_AM(u_ll, u_rr, orientation, equations)

		lambda = max(abs(v1_ll), abs(v1_rr))

		return f_ec - 0.5f0 * lambda * (u_rr - u_ll)
	end

	## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation
	@inline function flux_theta_rhos(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		# Unpack left and right state
		rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
		rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
		_, _, rho_theta_ll = u_ll
		_, _, rho_theta_rr = u_rr
		# Compute the necessary mean values
		rho_mean = ln_mean(rho_ll, rho_rr)
		gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)

		#gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

		v1_avg = 0.5f0 * (v1_ll + v1_rr)
		p_avg = 0.5f0 * (p_ll + p_rr)

		# Calculate fluxes depending on normal_direction
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr) * f1
		return SVector(f1, f2, f3)
	end

	@inline function flux_theta_rhos_es(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		_, v1_ll, _ = cons2prim(u_ll, equations)
		_, v1_rr, _ = cons2prim(u_rr, equations)

		f_ec = flux_theta_rhos(u_ll, u_rr, orientation, equations)

		lambda = max(abs(v1_ll), abs(v1_rr))

		return f_ec - 0.5f0 * lambda * (u_rr - u_ll)
	end

	@inline function flux_theta_rhos_AM(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		# Unpack left and right state
		rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
		rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
		_, _, rho_theta_ll = u_ll
		_, _, rho_theta_rr = u_rr
		# Compute the necessary mean values
		rho_mean = (rho_ll + rho_rr) * 0.5f0
		gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)

		#    gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

		v1_avg = 0.5f0 * (v1_ll + v1_rr)
		p_avg = 0.5f0 * (p_ll + p_rr)

		# Calculate fluxes depending on normal_direction
		f1 = rho_mean * v1_avg
		f2 = f1 * v1_avg + p_avg
		f3 = inv_ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr) * f1
		return SVector(f1, f2, f3)
	end

	@inline function flux_theta_rhos_AM_es(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		_, v1_ll, _ = cons2prim(u_ll, equations)
		_, v1_rr, _ = cons2prim(u_rr, equations)

		f_ec = flux_theta_rhos_AM(u_ll, u_rr, orientation, equations)

		lambda = max(abs(v1_ll), abs(v1_rr))

		return f_ec - 0.5f0 * lambda * (u_rr - u_ll)
	end

	## Entropy (total energy) conservative flux for the Compressible Euler with the Potential Formulation
	@inline function flux_theta_global(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		# Unpack left and right state
		rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
		rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)
		_, _, rho_theta_ll = u_ll
		_, _, rho_theta_rr = u_rr
		# Compute the necessary mean values
		rho_mean = ln_mean(rho_ll, rho_rr)
		gammamean = stolarsky_mean_opt(rho_theta_ll, rho_theta_rr, equations.gamma, p_ll, p_rr)

		#gammamean = stolarsky_mean(rho_theta_ll, rho_theta_rr, equations.gamma)

		v1_avg = 0.5f0 * (v1_ll + v1_rr)
		p_avg = 0.5f0 * (p_ll + p_rr)

		# Calculate fluxes depending on normal_direction
		f3 = gammamean * v1_avg
		f1 = f3 * ln_mean(rho_ll / rho_theta_ll, rho_rr / rho_theta_rr)
		f2 = f1 * v1_avg + p_avg
		return SVector(f1, f2, f3)
	end

	@inline function flux_theta_global_es(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		# Unpack left and right state
		_, v1_ll, _ = cons2prim(u_ll, equations)
		_, v1_rr, _ = cons2prim(u_rr, equations)

		f_ec = flux_theta_global(u_ll, u_rr, orientation, equations)

		lambda = max(abs(v1_ll), abs(v1_rr))

		return f_ec - 0.5f0 * lambda * (u_rr - u_ll)
	end

	@inline function prim2cons(prim, equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho, v1, p = prim
		rho_v1 = rho * v1
		rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
		return SVector(rho, rho_v1, rho_theta)
	end

	@inline function cons2prim(u, equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho, rho_v1, rho_theta = u
		v1 = rho_v1 / rho
		#  p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
		p = equations.K * (rho_theta)^equations.gamma
		p = equations.K * exp(log(rho_theta^equations.gamma))
		return SVector(rho, v1, p)
	end

	@inline function cons2cons(u, equations::CompressibleEulerPotentialTemperatureEquations1D)
		return u
	end

	@inline function cons2entropy2(u, equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho, rho_v1, rho_theta = u

		k = equations.gamma
		w1 = -0.5f0 * rho_v1^2 / (rho)^2
		w2 = rho_v1 / rho
		w3 = equations.gamma / (equations.gamma - 1) * k * (rho_theta)^(equations.gamma - 1)

		return SVector(w1, w2, w3)
	end

	@inline function cons2entropy(u, equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho, rho_v1, rho_theta = u

		k = equations.p_0 * (equations.R / equations.p_0)^equations.gamma
		w1 = log(k * (rho_theta / rho)^equations.gamma) - equations.gamma
		w2 = 0.0
		w3 = rho / rho_theta * equations.gamma

		return SVector(w1, w2, w3)
	end

	@inline function energy_total(cons, equations::CompressibleEulerPotentialTemperatureEquations1D)
		# Mathematical entropy
		p = equations.p_0 * (equations.R * cons[3] / equations.p_0)^equations.gamma

		U = (p / (equations.gamma - 1) + 1 / 2 * (cons[2]^2) / (cons[1]))

		return U
	end

	# Default entropy is the mathematical entropy
	@inline function entropy(cons, equations::CompressibleEulerPotentialTemperatureEquations1D)
		entropy_phys(cons, equations)
	end

	@inline function entropy_phys(cons, equations::CompressibleEulerPotentialTemperatureEquations1D)

		p = equations.p_0 * (equations.R * cons[3] / equations.p_0)^equations.gamma
		# Thermodynamic entropy
		s = log(p) - equations.gamma * log(cons[1])
		S = -s * cons[1] / (equations.gamma - 1.0)
		return S
	end

	@inline function energy_kinetic(cons, equations::CompressibleEulerPotentialTemperatureEquations1D)
		return 0.5f0 * (cons[2]^2) / (cons[1])
	end

	@inline function max_abs_speeds(u, equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho, v1, p = cons2prim(u, equations)
		c = sqrt(equations.gamma * p / rho)

		return (abs(v1) + c,)
	end

	@inline function max_abs_speed_naive(u_ll, u_rr, orientation::Integer,
		equations::CompressibleEulerPotentialTemperatureEquations1D)
		rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
		rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

		# Calculate primitive variables and speed of sound
		v_mag_ll = abs(v1_ll)
		c_ll = sqrt(equations.gamma * p_ll / rho_ll)
		v_mag_rr = abs(v1_rr)
		c_rr = sqrt(equations.gamma * p_rr / rho_rr)

		λ_max = max(v_mag_ll, v_mag_rr) + max(c_ll, c_rr)
	end

	@inline function pressure(cons, equations::CompressibleEulerPotentialTemperatureEquations1D)
		_, _, p = cons2prim(cons, equations)
		return p

	end


end # @muladd