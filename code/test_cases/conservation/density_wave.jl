using PotentialTemperature
using PotentialTemperature.Trixi
using DelimitedFiles
using CairoMakie
using LaTeXStrings
using OrdinaryDiffEqSSPRK

function setup_problem_density_wave(; polydeg::Int, dt::Float64, initial_refinement_level::Int64, surface_flux, volume_flux, problem_setup::ProblemSetup, use_volume_flux::Bool)

	@unpack tspan, periodicity, coordinates_min, coordinates_max, equations = problem_setup
	@unpack initial_condition, source_terms, boundary_conditions, time_method = problem_setup
	@unpack problem_name = problem_setup
	if use_volume_flux
		solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux, volume_integral = VolumeIntegralFluxDifferencing(volume_flux))
	else
		solver = DGSEM(polydeg = polydeg, surface_flux = surface_flux)
	end

	mesh = TreeMesh(coordinates_min, coordinates_max,
		initial_refinement_level = initial_refinement_level,
		n_cells_max = 100_000)
	semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver)

	ode = semidiscretize(semi, tspan)
	analysis_interval = Int(7500 / 4)
	summary_callback = SummaryCallback()
	if problem_name == "Euler"
		analysis_callback =
			AnalysisCallback(semi, interval = analysis_interval,
				save_analysis = true,
				output_directory = pwd() * "/test_cases/conservation/out/density_wave",
				analysis_filename = "analysiseuler.dat",
				extra_analysis_integrals = (energy_kinetic, energy_total, entropy, pressure))
	elseif problem_name == "Potential"
		analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
			save_analysis = true,
			output_directory = pwd() * "/test_cases/conservation/out/density_wave",
			analysis_filename = "analysistheta.dat",
			extra_analysis_integrals = (energy_kinetic, energy_total, entropy, pressure),
		)
	end
	alive_callback = AliveCallback(analysis_interval = 150000)

	stepsize_callback = StepsizeCallback(cfl = 0.01)

	callbacks = CallbackSet(summary_callback,
		analysis_callback,
		alive_callback,
		stepsize_callback)

	return ode, callbacks, semi

end

function initial_condition_density_wave(x, t, equations::CompressibleEulerPotentialTemperatureEquations1D)
	v1 = 1
	rho = 1 + exp(sinpi(2 * x[1]))
	rho_v1 = rho * v1
	p = 1
	rho_theta = (p / equations.p_0)^(1 / equations.gamma) * equations.p_0 / equations.R
	return SVector(rho, rho_v1, rho_theta)
end

function initial_condition_density_wave(x, t, equations::CompressibleEulerEquations1D)
	v1 = 1
	rho = 1 + exp(sinpi(2 * x[1]))
	rho_v1 = rho * v1
	p = 1
	rho_e = p / (equations.gamma - 1) + 0.5f0 * rho * v1^2
	return SVector(rho, rho_v1, rho_e)
end

function set_density_wave_1d(time_method, equations, problem_name)

	return ProblemSetup(problem_name = problem_name, equations = equations,
		initial_condition = initial_condition_density_wave,
		boundary_conditions = nothing,
		tspan = (0.0, 40.0),
		coordinates_min = (0.0,), coordinates_max = (1.0,),
		source_terms = nothing,
		periodicity = (true, true),
		time_method = time_method)
end

function run_density_wave_1d(; time_method = SSPRK43(), polydeg::Int = 0, initial_refinement_level = 6)

	EulerTGV = set_density_wave_1d(time_method, CompressibleEulerEquations1D(1004.0 / 717.0), "Euler")
	PotentialTGV = set_density_wave_1d(time_method, CompressibleEulerPotentialTemperatureEquations1D(), "Potential")

	fluxes = (flux_theta, flux_theta_AM, flux_theta_rhos, flux_theta_rhos_AM, flux_theta_global)
	density_wave_1d(; polydeg = polydeg, time_method = time_method, initial_refinement_level = initial_refinement_level, surface_flux = flux_ranocha, volume_flux = flux_ranocha, problem_setup = EulerTGV, use_volume_flux = true)
	for flux in fluxes
		density_wave_1d(; polydeg = polydeg, time_method = time_method, initial_refinement_level = initial_refinement_level, surface_flux = flux, volume_flux = flux, problem_setup = PotentialTGV, use_volume_flux = true)
	end
end

function density_wave_1d(; polydeg::Int, dt = 1.0, time_method, initial_refinement_level::Int64, surface_flux, volume_flux, problem_setup::ProblemSetup, use_volume_flux::Bool)
	@unpack equations, problem_name = problem_setup

	ode, callbacks, semi = setup_problem_density_wave(; polydeg = polydeg, dt = dt, initial_refinement_level = initial_refinement_level,
		surface_flux = surface_flux, volume_flux = volume_flux,
		problem_setup = problem_setup, use_volume_flux = use_volume_flux)

	try
		sol = solve(ode, time_method,
			dt = 1.0, # solve needs some value here but it will be overwritten by the stepsize_callback
			save_everystep = false, callback = callbacks, adaptive = false)
	catch error
		@info "Blow-up" error
	end


	if problem_name == "Euler"
		data = readdlm(pwd() * "/test_cases/conservation/out/density_wave/analysiseuler.dat", skipstart = 1)
	elseif problem_name == "Potential"
		data = readdlm(pwd() * "/test_cases/conservation/out/density_wave/analysistheta.dat", skipstart = 1)
	end

	col_time = 2                 # Time column
	col_energy_kinetic = 11      # Energy Kinetic
	col_total_energy = 12        # Entropy
	col_entropy_phys = 13        # Entropy Phys
	col_pressure = 14            # Pressure

	time = data[:, col_time]
	energy_kinetic = (data[:, col_energy_kinetic] .- data[1, col_energy_kinetic]) ./ data[1, col_energy_kinetic]
	total_energy = (data[:, col_total_energy] .- data[1, col_total_energy]) ./ data[1, col_total_energy]
	entropy_phys = (data[:, col_entropy_phys] .- data[1, col_entropy_phys]) ./ data[1, col_entropy_phys]
	pressure = (data[:, col_pressure] .- data[1, col_pressure]) ./ data[1, col_pressure]

	results = hcat(time, energy_kinetic, total_energy, entropy_phys, pressure)

	if problem_name == "Euler"
		writedlm(pwd() * "/test_cases/conservation/out/density_wave/euler_ref$(initial_refinement_level)_polydeg$polydeg.dat", results)
	elseif problem_name == "Potential"
		writedlm(pwd() * "/test_cases/conservation/out/density_wave/potential_" * string(surface_flux) * "_ref$(initial_refinement_level)_polydeg$polydeg.dat", results)
	end

	return nothing
end

run_density_wave_1d()