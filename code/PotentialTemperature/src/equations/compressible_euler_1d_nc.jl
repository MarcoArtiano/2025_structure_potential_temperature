using Trixi
using Trixi: ln_mean, stolarsky_mean, AbstractCompressibleEulerEquations
import Trixi: varnames, cons2cons, cons2prim, cons2entropy, entropy, have_nonconservative_terms, boundary_condition_slip_wall, flux_ranocha, flux_shima_etal, flux_kennedy_gruber, flux_chandrashekar, flux_hllc
# By default, Julia/LLVM does not use fused multiply-add operations (FMAs).
# Since these FMAs can increase the performance of many numerical algorithms,
# we need to opt-in explicitly.
# See https://ranocha.de/blog/Optimizing_EC_Trixi for further details.
@muladd begin
#! format: noindent

@doc raw"""
	CompressibleEulerEquations1D(gamma)

The compressible Euler equations
```math
\frac{\partial}{\partial t}
\begin{pmatrix}
\rho \\ \rho v_1 \\ \rho e
\end{pmatrix}
+
\frac{\partial}{\partial x}
\begin{pmatrix}
\rho v_1 \\ \rho v_1^2 + p \\ (\rho e +p) v_1
\end{pmatrix}
=
\begin{pmatrix}
0 \\ 0 \\ 0
\end{pmatrix}
```
for an ideal gas with ratio of specific heats `gamma` in one space dimension.
Here, ``\rho`` is the density, ``v_1`` the velocity, ``e`` the specific total energy **rather than** specific internal energy, and
```math
p = (\gamma - 1) \left( \rho e - \frac{1}{2} \rho v_1^2 \right)
```
the pressure.
"""
struct CompressibleEulerEquations1DNC{RealT <: Real} <:
	   AbstractCompressibleEulerEquations{1, 4}
	gamma::RealT               # ratio of specific heats
	inv_gamma_minus_one::RealT # = inv(gamma - 1); can be used to write slow divisions as fast multiplications

	function CompressibleEulerEquations1DNC(gamma)
		γ, inv_gamma_minus_one = promote(gamma, inv(gamma - 1))
		new{typeof(γ)}(γ, inv_gamma_minus_one)
	end
end

function varnames(::typeof(cons2cons), ::CompressibleEulerEquations1DNC)
	("rho", "rho_v1", "rho_e", "phi")
end
varnames(::typeof(cons2prim), ::CompressibleEulerEquations1DNC) = ("rho", "v1", "p", "phi")

have_nonconservative_terms(::CompressibleEulerEquations1DNC) = Trixi.True()


@inline function boundary_condition_slip_wall_2(u_inner, orientation_or_normal, direction,
	x, t,
	surface_flux_function,
	equations::CompressibleEulerEquations1DNC)

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

# Calculate 1D flux for a single point
@inline function flux(u, orientation::Integer, equations::CompressibleEulerEquations1DNC)
	rho, rho_v1, rho_e, phi = u
	v1 = rho_v1 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho_v1 * v1)
	# Ignore orientation since it is always "1" in 1D
	f1 = rho_v1
	f2 = rho_v1 * v1 + p
	f3 = (rho_e + p) * v1
	return SVector(f1, f2, f3, 0)
end

@inline function flux_kennedy_gruber(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations1DNC)
	# Unpack left and right state
	_, _, rho_e_ll, _ = u_rr
	_, _, rho_e_rr, _ = u_rr
	rho_ll, v1_ll, p_ll = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr = cons2prim(u_rr, equations)

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	e_avg = 0.5f0 * (rho_e_ll / rho_ll + rho_e_rr / rho_rr)

	# Ignore orientation since it is always "1" in 1D
	f1 = rho_avg * v1_avg
	f2 = rho_avg * v1_avg * v1_avg + p_avg
	f3 = (rho_avg * e_avg + p_avg) * v1_avg

	return SVector(f1, f2, f3, 0.0)
end

"""
	flux_chandrashekar(u_ll, u_rr, orientation, equations::CompressibleEulerEquations1D)

Entropy conserving two-point flux by
- Chandrashekar (2013)
  Kinetic Energy Preserving and Entropy Stable Finite Volume Schemes
  for Compressible Euler and Navier-Stokes Equations
  [DOI: 10.4208/cicp.170712.010313a](https://doi.org/10.4208/cicp.170712.010313a)
"""
@inline function flux_chandrashekar(u_ll, u_rr, orientation::Integer,
	equations::CompressibleEulerEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, phi = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, phi = cons2prim(u_rr, equations)
	beta_ll = 0.5f0 * rho_ll / p_ll
	beta_rr = 0.5f0 * rho_rr / p_rr
	specific_kin_ll = 0.5f0 * (v1_ll^2)
	specific_kin_rr = 0.5f0 * (v1_rr^2)

	# Compute the necessary mean values
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	rho_mean = ln_mean(rho_ll, rho_rr)
	beta_mean = ln_mean(beta_ll, beta_rr)
	beta_avg = 0.5f0 * (beta_ll + beta_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_mean = 0.5f0 * rho_avg / beta_avg
	velocity_square_avg = specific_kin_ll + specific_kin_rr

	# Calculate fluxes
	# Ignore orientation since it is always "1" in 1D
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_mean
	f3 = f1 * 0.5f0 * (1 / (equations.gamma - 1) / beta_mean - velocity_square_avg) +
		 f2 * v1_avg

	return SVector(f1, f2, f3, 0.0)
end


"""
	flux_shima_etal(u_ll, u_rr, orientation, equations::CompressibleEulerEquations1D)

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
	equations::CompressibleEulerEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, phi = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, phi = cons2prim(u_rr, equations)

	# Average each factor of products in flux
	rho_avg = 0.5f0 * (rho_ll + rho_rr)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	kin_avg = 0.5f0 * (v1_ll * v1_rr)

	# Calculate fluxes
	# Ignore orientation since it is always "1" in 1D
	pv1_avg = 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)
	f1 = rho_avg * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = p_avg * v1_avg * equations.inv_gamma_minus_one + f1 * kin_avg + pv1_avg

	return SVector(f1, f2, f3, 0)
end


"""
	flux_ranocha(u_ll, u_rr, orientation_or_normal_direction, equations::CompressibleEulerEquations1D)

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
	equations::CompressibleEulerEquations1DNC)
	# Unpack left and right state
	rho_ll, v1_ll, p_ll, phi = cons2prim(u_ll, equations)
	rho_rr, v1_rr, p_rr, phi = cons2prim(u_rr, equations)

	# Compute the necessary mean values
	rho_mean = ln_mean(rho_ll, rho_rr)
	# Algebraically equivalent to `inv_ln_mean(rho_ll / p_ll, rho_rr / p_rr)`
	# in exact arithmetic since
	#     log((ϱₗ/pₗ) / (ϱᵣ/pᵣ)) / (ϱₗ/pₗ - ϱᵣ/pᵣ)
	#   = pₗ pᵣ log((ϱₗ pᵣ) / (ϱᵣ pₗ)) / (ϱₗ pᵣ - ϱᵣ pₗ)
	inv_rho_p_mean = p_ll * p_rr * inv_ln_mean(rho_ll * p_rr, rho_rr * p_ll)
	v1_avg = 0.5f0 * (v1_ll + v1_rr)
	p_avg = 0.5f0 * (p_ll + p_rr)
	velocity_square_avg = 0.5f0 * (v1_ll * v1_rr)

	# Calculate fluxes
	# Ignore orientation since it is always "1" in 1D
	f1 = rho_mean * v1_avg
	f2 = f1 * v1_avg + p_avg
	f3 = f1 * (velocity_square_avg + inv_rho_p_mean * equations.inv_gamma_minus_one) +
		 0.5f0 * (p_ll * v1_rr + p_rr * v1_ll)

	return SVector(f1, f2, f3, 0.0)
end

# While `normal_direction` isn't strictly necessary in 1D, certain solvers assume that
# the normal component is incorporated into the numerical flux.
#
# See `flux(u, normal_direction::AbstractVector, equations::AbstractEquations{1})` for a
# similar implementation.
@inline function flux_ranocha_nc(u_ll, u_rr, normal_direction::AbstractVector,
	equations::CompressibleEulerEquations1DNC)
	return normal_direction[1] * flux_ranocha(u_ll, u_rr, 1, equations)
end

@inline function max_abs_speeds(u, equations::CompressibleEulerEquations1DNC)
	rho, rho_v1, rho_e, phi = u
	v1 = rho_v1 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * v1^2)
	c = sqrt(equations.gamma * p / rho)

	return (abs(v1) + c,)
end

# Convert conservative variables to primitive
@inline function cons2prim(u, equations::CompressibleEulerEquations1DNC)
	rho, rho_v1, rho_e, phi = u

	v1 = rho_v1 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho_v1 * v1)

	return SVector(rho, v1, p, phi)
end

# Convert conservative variables to entropy
@inline function cons2entropy(u, equations::CompressibleEulerEquations1DNC)
	rho, rho_v1, rho_e, phi = u

	v1 = rho_v1 / rho
	v_square = v1^2
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * rho * v_square)
	s = log(p) - equations.gamma * log(rho)
	rho_p = rho / p

	w1 = (equations.gamma - s) * equations.inv_gamma_minus_one -
		 0.5f0 * rho_p * v_square
	w2 = rho_p * v1
	w3 = -rho_p

	return SVector(w1, w2, w3, phi)
end

@inline function entropy2cons(w, equations::CompressibleEulerEquations1DNC)
	# See Hughes, Franca, Mallet (1986) A new finite element formulation for CFD
	# [DOI: 10.1016/0045-7825(86)90127-1](https://doi.org/10.1016/0045-7825(86)90127-1)
	@unpack gamma = equations

	# convert to entropy `-rho * s` used by Hughes, France, Mallet (1986)
	# instead of `-rho * s / (gamma - 1)`
	V1, V2, V5, phi = w .* (gamma - 1)

	# specific entropy, eq. (53)
	s = gamma - V1 + 0.5f0 * (V2^2) / V5

	# eq. (52)
	energy_internal = ((gamma - 1) / (-V5)^gamma)^(equations.inv_gamma_minus_one) *
					  exp(-s * equations.inv_gamma_minus_one)

	# eq. (51)
	rho = -V5 * energy_internal
	rho_v1 = V2 * energy_internal
	rho_e = (1 - 0.5f0 * (V2^2) / V5) * energy_internal
	return SVector(rho, rho_v1, rho_e, phi)
end

# Convert primitive to conservative variables
@inline function prim2cons(prim, equations::CompressibleEulerEquations1DNC)
	rho, v1, p, phi = prim
	rho_v1 = rho * v1
	rho_e = p * equations.inv_gamma_minus_one + 0.5f0 * (rho_v1 * v1)
	return SVector(rho, rho_v1, rho_e, phi)
end

@inline function density(u, equations::CompressibleEulerEquations1DNC)
	rho = u[1]
	return rho
end

@inline function velocity(u, equations::CompressibleEulerEquations1DNC)
	rho = u[1]
	v1 = u[2] / rho
	return v1
end

@inline function pressure(u, equations::CompressibleEulerEquations1DNC)
	rho, rho_v1, rho_e, phi = u
	p = (equations.gamma - 1) * (rho_e - 0.5f0 * (rho_v1^2) / rho)
	return p
end

@inline function density_pressure(u, equations::CompressibleEulerEquations1DNC)
	rho, rho_v1, rho_e, phi = u
	rho_times_p = (equations.gamma - 1) * (rho * rho_e - 0.5f0 * (rho_v1^2))
	return rho_times_p
end

# Calculate thermodynamic entropy for a conservative state `cons`
@inline function entropy_thermodynamic(cons, equations::CompressibleEulerEquations1D)
	# Pressure
	p = (equations.gamma - 1) * (cons[3] - 0.5f0 * (cons[2]^2) / cons[1])

	# Thermodynamic entropy
	s = log(p) - equations.gamma * log(cons[1])

	return s
end

# Calculate mathematical entropy for a conservative state `cons`
@inline function entropy_math(cons, equations::CompressibleEulerEquations1DNC)
	# Mathematical entropy
	S = -entropy_thermodynamic(cons, equations) * cons[1] *
		equations.inv_gamma_minus_one

	return S
end

# Default entropy is the mathematical entropy
@inline function entropy(cons, equations::CompressibleEulerEquations1DNC)
	entropy_math(cons, equations)
end

# Calculate total energy for a conservative state `cons`
@inline energy_total(cons, ::CompressibleEulerEquations1DNC) = cons[3]

# Calculate kinetic energy for a conservative state `cons`
@inline function energy_kinetic(cons, equations::CompressibleEulerEquations1DNC)
	return 0.5f0 * (cons[2]^2) / cons[1]
end

# Calculate internal energy for a conservative state `cons`
@inline function energy_internal(cons, equations::CompressibleEulerEquations1DNC)
	return energy_total(cons, equations) - energy_kinetic(cons, equations)
end
end # @muladd
