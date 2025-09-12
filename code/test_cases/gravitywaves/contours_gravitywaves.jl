using OrdinaryDiffEqSSPRK
using PotentialTemperature.Trixi
using PotentialTemperature
using CairoMakie

function initial_condition_gravity_wave(x, t,
	equations::CompressibleEulerPotentialTemperatureEquations2DNC)
	g = equations.g
	c_p = equations.c_p
	c_v = equations.c_v
	# center of perturbation
	x_c = 100_000.0
	a = 5_000
	H = 10_000
	# perturbation in potential temperature
	R = c_p - c_v    # gas constant (dry air)

	T0 = 250
	delta = 9.81 / (R * T0)
	DeltaT = 0.001
	Tb = DeltaT * sinpi(x[2] / H) * exp(-(x[1] - x_c)^2 / a^2)
	ps = 100_000.0  # reference pressure
	rhos = ps / (T0 * R)
	rho_b = rhos * (-Tb / T0)
	p = ps * exp(-delta * x[2])
	rho = rhos * exp(-delta * x[2]) + rho_b * exp(-0.5 * delta * x[2])
	v1 = 20.0
	v2 = 0.0

	return prim2cons(SVector(rho, v1, v2, p, x[2]), equations)
end

function initial_condition_gravity_wave(x, t,
	equations::CompressibleEulerPotentialTemperatureEquations2D)
	g = equations.g
	c_p = equations.c_p
	c_v = equations.c_v
	# center of perturbation
	x_c = 100_000.0
	a = 5_000
	H = 10_000
	# perturbation in potential temperature
	R = c_p - c_v    # gas constant (dry air)

	T0 = 250
	delta = 9.81 / (R * T0)
	DeltaT = 0.001
	Tb = DeltaT * sinpi(x[2] / H) * exp(-(x[1] - x_c)^2 / a^2)
	ps = 100_000.0  # reference pressure
	rhos = ps / (T0 * R)
	rho_b = rhos * (-Tb / T0)
	p = ps * exp(-delta * x[2])
	rho = rhos * exp(-delta * x[2]) + rho_b * exp(-0.5 * delta * x[2])
	v1 = 20.0
	v2 = 0.0

	return prim2cons(SVector(rho, v1, v2, p), equations)
end
@inline function flux_zero(u_ll, u_rr, normal_direction, equations)
    return zero(u_ll)
end
equations = CompressibleEulerPotentialTemperatureEquations2DNC()
cs = sqrt(1004.0 / 717.0 * 287.0 * 250.0)
surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_theta, flux_nonconservative_gravity_log)
polydeg = 3
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux, volume_integral = VolumeIntegralFluxDifferencing(volume_flux))

boundary_conditions = (x_neg = boundary_condition_periodic,
	x_pos = boundary_condition_periodic,
	y_neg = boundary_condition_slip_wall_2,
	y_pos = boundary_condition_slip_wall_2)

coordinates_min = (0.0, 0.0)
coordinates_max = (300_000.0, 10_000.0)
cells_per_dimension = (30*3, 4*3)
mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max, periodicity = (true, false))
source_terms = nothing
initial_condition = initial_condition_gravity_wave
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms = source_terms,
	boundary_conditions = boundary_conditions)
tspan = (0.0, 1800.0)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
	analysis_callback,
	alive_callback,
	stepsize_callback)

time_method = SSPRK43()
sol = solve(ode,
	time_method,
	maxiters = 1.0e7,
	dt = 1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
	save_everystep = false, callback = callbacks, adaptive = false)

x, y, data = ContourData(sol.u[end], semi, cells_per_dimension, equations)

h = Figure(size = (800, 300))
labelsize = 20
kwargs = (xlabel = L"$x$ [km]", xlabelsize = labelsize, ylabelsize = labelsize, limits = ((0, 300), (0, 10)), xticklabelsize = 17.0, yticklabelsize = 17.0, ylabel = L"$z$ [km]")

Axis(h[1, 1]; kwargs...)

c = contourf!(x ./ 1e3, y ./ 1e3, data[3, :, :] ./ data[1, :, :], levels = 7, colormap = :cividis)
Colorbar(h[1, 2], c, ticklabelsize = 17)

equations = CompressibleEulerPotentialTemperatureEquations2D()
surface_flux = FluxLMARS(cs)
volume_flux = flux_theta
solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux, volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
source_terms = source_terms_gravity
semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms = source_terms,
	boundary_conditions = boundary_conditions)

ode = semidiscretize(semi, tspan)
analysis_callback = AnalysisCallback(semi, interval = analysis_interval)

alive_callback = AliveCallback(analysis_interval = analysis_interval)

stepsize_callback = StepsizeCallback(cfl = 1.0)

callbacks = CallbackSet(summary_callback,
	analysis_callback,
	alive_callback,
	stepsize_callback)
sol = solve(ode,
	time_method,
	maxiters = 1.0e7,
	dt = 1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
	save_everystep = false, callback = callbacks, adaptive = false)

x, y, data = ContourData(sol.u[end], semi, cells_per_dimension, equations)

kwargs = (xlabel = L"$x$ [km]", xlabelsize = labelsize, ylabelsize = labelsize, limits = ((0, 300), (0, 10)), xticklabelsize = 17.0, yticklabelsize = 17.0, ylabel = L"$z$ [km]")

Axis(h[2, 1]; kwargs...)

c = contourf!(x ./ 1e3, y ./ 1e3, data[3, :, :] ./ data[1, :, :], levels = 7, colormap = :cividis)
Colorbar(h[2, 2], c, ticklabelsize = 17)

save(pwd() * "/test_cases/gravitywaves/plots/contour_comparison_$(cells_per_dimension[1])x$(cells_per_dimension[2])_CFL1_polydeg$(polydeg)_pointwise.pdf", h)