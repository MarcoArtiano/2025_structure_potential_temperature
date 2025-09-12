using OrdinaryDiffEqSSPRK
using PotentialTemperature.Trixi
using PotentialTemperature
using CairoMakie

struct HydrostaticSetup
	# Physical constants
	g::Float64       # gravity of earth
	c_p::Float64     # heat capacity for constant pressure (dry air)
	c_v::Float64     # heat capacity for constant volume (dry air)
	gamma::Float64   # heat capacity ratio (dry air)
	p_0::Float64     # atmospheric pressure
	T_0::Float64     # 
	u0::Float64      #
	Nf::Float64      # 
	z_B::Float64     #
	z_T::Float64     #
	alfa::Float64    #
	xr_B::Float64
	form1::Bool
	form2::Bool
	function HydrostaticSetup(alfa, xr_B, form1, form2; g = 9.81, c_p = 1004.0, c_v = 717.0, gamma = c_p / c_v, p_0 = 100_000.0, T_0 = 250.0, u0 = 20.0, z_B = 15000.0, z_T = 30000.0)
		Nf = g / sqrt(c_p * T_0)
		new(g, c_p, c_v, gamma, p_0, T_0, u0, Nf, z_B, z_T, alfa, xr_B, form1, form2)
	end
end


@inline function rayleigh_damping(x, z_B, z_T, alfa, xr_B)
	xr_T = 120000.0

	if x[2] <= z_B
		S_v = 0.0
	else
		S_v = -alfa * sinpi(0.5 * (x[2] - z_B) / (z_T - z_B))^2
	end
	if x[1] < xr_B
		S_h1 = 0.0
	else
		S_h1 = -alfa * sinpi(0.5 * (x[1] - xr_B) / (xr_T - xr_B))^2
	end

	if x[1] > -xr_B
		S_h2 = 0.0
	else
		S_h2 = -alfa * sinpi(0.5 * (x[1] + xr_B) / (-xr_T + xr_B))^2
	end
	return S_v, S_h1, S_h2
end

function (setup::HydrostaticSetup)(u, x, t, equations::CompressibleEulerEquations2DNC)
	@unpack g, c_p, c_v, gamma, p_0, T_0, z_B, z_T, Nf, u0, alfa, xr_B, form1, form2 = setup

	rho, rho_v1, rho_v2, rho_e, _ = u

	R = c_p - c_v

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	rho_theta = (p / p_0)^(c_v / c_p) * p_0 / R
	theta = rho_theta / rho

	S_v, S_h1, S_h2 = rayleigh_damping(x, z_B, z_T, alfa, xr_B)

	exner = exp(-Nf^2 / g * x[2])

	theta_0 = T_0 / exner
	p00 = p_0 * exner^(c_p / R)
	rho00 = p00 / (R * T_0)
	rho_E_0 = p00 / (equations.gamma - 1) + 0.5 * rho00 * 20^2

	K = p_0 * (R / p_0)^gamma
	du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
	du3 = rho_v2 * (S_v + S_h1 + S_h2)
	du42 = rho * (theta - theta_0) * (S_v + S_h1 + S_h2) * K * gamma / (gamma - 1.0) * (rho_theta)^(gamma - 1.0) + du2 * v1 + du3 * v2
	du41 = (rho_e - rho_E_0) * (S_v + S_h1 + S_h2)
	return SVector(zero(eltype(u)), du2, du3, du41*form1 + du42*form2, zero(eltype(u)))

end

function (setup::HydrostaticSetup)(u, x, t, equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	@unpack g, c_p, c_v, gamma, p_0, T_0, z_B, z_T, Nf, u0, alfa, xr_B = setup

	rho, rho_v1, rho_v2, rho_theta, _ = u

	R = c_p - c_v

	v1 = rho_v1 / rho
	theta = rho_theta / rho

	S_v, S_h1, S_h2 = rayleigh_damping(x, z_B, z_T, alfa, xr_B)

	exner = exp(-Nf^2 / g * x[2])

	theta_0 = T_0 / exner

	du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
	du3 = rho_v2 * (S_v + S_h1 + S_h2)
	du4 = rho * (theta - theta_0) * (S_v + S_h1 + S_h2)

	return SVector(zero(eltype(u)), du2, du3, du4, zero(eltype(u)))

end

function (setup::HydrostaticSetup)(u, x, t, equations::CompressibleEulerEquations2D)
	@unpack g, c_p, c_v, gamma, p_0, T_0, z_B, z_T, Nf, u0, alfa, xr_B, form1, form2 = setup

	rho, rho_v1, rho_v2, rho_e = u

	R = c_p - c_v

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	rho_theta = (p / p_0)^(c_v / c_p) * p_0 / R
	theta = rho_theta / rho

	S_v, S_h1, S_h2 = rayleigh_damping(x, z_B, z_T, alfa, xr_B)

	exner = exp(-Nf^2 / g * x[2])

	theta_0 = T_0 / exner

	p_0 = 100_000.0  # reference pressure
	R = c_p - c_v    # gas constant (dry air)
	p00 = p_0 * exner^(c_p / R)

	# density
	rho00 = p00 / (R * T_0)

	rho_E_0 = p00 / (equations.gamma - 1) + 0.5 * rho00 * 20^2
	K = p_0 * (R / p_0)^gamma
	du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
	du3 = rho_v2 * (S_v + S_h1 + S_h2)
	du42 = rho * (theta - theta_0) * (S_v + S_h1 + S_h2) * K * gamma / (gamma - 1.0) * (rho_theta)^(gamma - 1.0) + du2 * v1 + du3 * v2
	du41 = (rho_e - rho_E_0) * (S_v + S_h1 + S_h2)
	
	return SVector(zero(eltype(u)), du2, du3 - g * rho, du41*form1 + du42*form2 - g * rho_v2)

end

function (setup::HydrostaticSetup)(u, x, t, equations::CompressibleEulerPotentialTemperatureEquations2D)
	@unpack g, c_p, c_v, gamma, p_0, T_0, z_B, z_T, Nf, u0, alfa, xr_B = setup

	rho, rho_v1, rho_v2, rho_theta = u

	R = c_p - c_v

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	theta = rho_theta / rho

	S_v, S_h1, S_h2 = rayleigh_damping(x, z_B, z_T, alfa, xr_B)

	exner = exp(-Nf^2 / g * x[2])

	theta_0 = T_0 / exner

	du2 = rho * (v1 - u0) * (S_v + S_h1 + S_h2)
	du3 = rho_v2 * (S_v + S_h1 + S_h2)
	du4 = rho * (theta - theta_0) * (S_v + S_h1 + S_h2)

	return SVector(zero(eltype(u)), du2, du3 - g * rho, du4)

end

function (setup::HydrostaticSetup)(x, t, equations::Union{CompressibleEulerEquations2D, CompressibleEulerPotentialTemperatureEquations2D})
	@unpack g, c_p, c_v, p_0, T_0, u0, Nf = setup

	# Exner pressure, solves hydrostatic equation for x[2]
	exner = exp(-Nf^2 / g * x[2])
	# pressure
	p_0 = 100_000.0  # reference pressure
	R = c_p - c_v    # gas constant (dry air)
	p = p_0 * exner^(c_p / R)

	# density
	rho = p / (R * T_0)
	v1 = u0
	v2 = 0.0

	return prim2cons(SVector(rho, v1, v2, p), equations)
end

function (setup::HydrostaticSetup)(x, t, equations::Union{CompressibleEulerEquations2DNC, CompressibleEulerPotentialTemperatureEquations2DNC})
	@unpack g, c_p, c_v, p_0, T_0, u0, Nf = setup

	# Exner pressure, solves hydrostatic equation for x[2]
	exner = exp(-Nf^2 / g * x[2])
	# pressure
	p_0 = 100_000.0  # reference pressure
	R = c_p - c_v    # gas constant (dry air)
	p = p_0 * exner^(c_p / R)

	# density
	rho = p / (R * T_0)
	v1 = u0
	v2 = 0.0

	return prim2cons(SVector(rho, v1, v2, p, x[2]), equations)
end

    equations = CompressibleEulerEquations2DNC(1007/717)
    alfa = 0.035
    form1 = false
    xr_B = 60000.0
	linear_hydrostatic_setup = HydrostaticSetup(alfa, xr_B, form1, !form1)

	boundary = BoundaryConditionDirichlet(linear_hydrostatic_setup)

	boundary_conditions = Dict(:x_neg => boundary,
		:x_pos => boundary,
		:y_neg => boundary_condition_slip_wall_2,
		:y_pos => boundary)
	cs = sqrt(1004.0 / 717.0 * 287.0 * 250.0)
    polydeg = 3
	basis = LobattoLegendreBasis(polydeg)
    surface_flux = (FluxLMARS(cs), flux_nonconservative_gravity_log)
    volume_flux = (flux_ranocha, flux_nonconservative_gravity_log)
	volume_integral = VolumeIntegralFluxDifferencing(volume_flux)

	solver = DGSEM(basis, surface_flux, volume_integral)

	a = 10000.0
	L = 240000.0
	H = 30000.0
	peak = 1.0
	y_b = peak / (1 + (L / 2 / a)^2)
	alfa_b = (H - y_b) * 0.5

	f1(s) = SVector(-L / 2, y_b + alfa_b * (s + 1))
	f2(s) = SVector(L / 2, y_b + alfa_b * (s + 1))
	f3(s) = SVector((s + 1 - 1) * L / 2, peak / (1 + ((s + 1 - 1) * L / 2)^2 / a^2))
	f4(s) = SVector((s + 1 - 1) * L / 2, H)
    cells_per_dimension = (100, 60)
    mesh = P4estMesh(cells_per_dimension, polydeg = polydeg,
		faces = (f1, f2, f3, f4),
		initial_refinement_level = 0, periodicity = (false, false))

	semi = SemidiscretizationHyperbolic(mesh, equations, linear_hydrostatic_setup, solver, source_terms = linear_hydrostatic_setup,
		boundary_conditions = boundary_conditions)

	###############################################################################
	# ODE solvers, callbacks etc.
    T = 12.5
	tspan = (0.0, T * 3600.0)  # 1000 seconds final time
	ode = semidiscretize(semi, tspan)

	summary_callback = SummaryCallback()

	analysis_interval = 1000

	analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

	alive_callback = AliveCallback(analysis_interval = analysis_interval)

	callbacks = CallbackSet(summary_callback,
		analysis_callback,
		alive_callback)

	###############################################################################
	# run the simulation
	sol = solve(ode,
		SSPRK43(),
		maxiters = 1.0e7,
		dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
		save_everystep = false, callback = callbacks)

x, y, data = ContourData(sol.u[end], semi, cells_per_dimension, equations)

h = Figure(size = (850, 400))
labelsize = 20
levels = 13
kwargs = (xlabel = L"$x$ [km]", xlabelsize = labelsize, ylabelsize = labelsize, limits = ((-40, 40), (0, 12)), xticklabelsize = 17.0, yticklabelsize = 17.0)

Axis(h[1, 1]; kwargs..., ylabel = L"$z$ [km]")

c = contourf!(x ./ 1e3, y ./ 1e3, data[2, :, :] ./ data[1, :, :], levels = levels, colormap = :cividis)
Colorbar(h[1, 2], c, ticklabelsize = 17)

Axis(h[1, 3]; kwargs...)

c = contourf!(x ./ 1e3, y ./ 1e3, data[3, :, :] ./ data[1, :, :], levels = levels, colormap = :cividis)
Colorbar(h[1, 4], c, ticklabelsize = 17)

save(pwd() * "/test_cases/linearhydrostatic/plots/contour_comparison_$(cells_per_dimension[1])x$(cells_per_dimension[2])_CFL1_polydeg$(polydeg)_u_w.pdf", h)