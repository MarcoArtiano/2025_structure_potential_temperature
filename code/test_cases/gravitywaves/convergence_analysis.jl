using OrdinaryDiffEqSSPRK
using PotentialTemperature.Trixi
using DataFrames
using PotentialTemperature
using CSV

include("inertia_gravity_analytical_solution.jl")

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

function initial_condition_gravity_wave(x, t,
    equations::CompressibleEulerEquations2DNC)
    g = 9.81
    c_p = 1004.0
    c_v = 717.0
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
    equations::CompressibleEulerEquations2D)
    g = 9.81
    c_p = 1004.0
    c_v = 717.0
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

function setup_problem_inertia_gravity_wave(; polydeg::Int, CFL::Float64, cells_per_dimension::Tuple{Int64,Int64}, surface_flux, volume_flux, problem_setup::ProblemSetup)

    @unpack tspan, periodicity, coordinates_min, coordinates_max, equations = problem_setup
    @unpack initial_condition, source_terms, boundary_conditions, time_method = problem_setup

    solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux, volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

    mesh = StructuredMesh(cells_per_dimension, coordinates_min, coordinates_max, periodicity=periodicity)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms=source_terms,
        boundary_conditions=boundary_conditions)

    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 10000
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    stepsize_callback = StepsizeCallback(cfl=CFL)

    callbacks = CallbackSet(summary_callback,
        analysis_callback,
        alive_callback,
        stepsize_callback)

    return ode, callbacks, semi

end

function setup_problem_inertia_gravity_wave_warped(; polydeg::Int, CFL::Float64, cells_per_dimension::Tuple{Int64,Int64}, surface_flux, volume_flux, problem_setup::ProblemSetup)

    @unpack tspan, periodicity, coordinates_min, coordinates_max, equations = problem_setup
    @unpack initial_condition, source_terms, boundary_conditions, time_method = problem_setup

    solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux, volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

    function mapping(xi, eta)

        H = 10000
        L = 300000
        xit = L * 0.5 * (1 + xi)
        etat = H * 0.5 * (1 + eta)  
        x = xit + L/20 * sinpi(xit/L) * sinpi(2* etat/H)
        y = etat - H/20 * sinpi(2*xit/L) * sinpi(  etat/H)
        return SVector(x, y)
    end

    mesh = StructuredMesh(cells_per_dimension, mapping, periodicity=periodicity)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver, source_terms=source_terms,
        boundary_conditions=boundary_conditions)

    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 10000
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    stepsize_callback = StepsizeCallback(cfl=CFL)

    callbacks = CallbackSet(summary_callback,
        analysis_callback,
        alive_callback,
        stepsize_callback)

    return ode, callbacks, semi

end

function inertia_gravity_wave(; polydeg::Int, CFL::Float64, time_method, cells_per_dimension::Tuple{Int64,Int64}, surface_flux, volume_flux, problem_setup::ProblemSetup)
    @unpack equations = problem_setup

    ode, callbacks, semi = setup_problem_inertia_gravity_wave(; polydeg=polydeg, CFL=CFL, cells_per_dimension=cells_per_dimension,
        surface_flux=surface_flux, volume_flux=volume_flux,
        problem_setup=problem_setup)

    sol = solve(ode,
        time_method,
        maxiters=1.0e7,
        dt=1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
        save_everystep=false, callback=callbacks, adaptive=false)

    return sol, semi

end

function inertia_gravity_wave_warped(; polydeg::Int, CFL::Float64, time_method, cells_per_dimension::Tuple{Int64,Int64}, surface_flux, volume_flux, problem_setup::ProblemSetup)
    @unpack equations = problem_setup

    ode, callbacks, semi = setup_problem_inertia_gravity_wave_warped(; polydeg=polydeg, CFL=CFL, cells_per_dimension=cells_per_dimension,
        surface_flux=surface_flux, volume_flux=volume_flux,
        problem_setup=problem_setup)

    sol = solve(ode,
        time_method,
        maxiters=1.0e7,
        dt=1e-1, # solve needs some value here but it will be overwritten by the stepsize_callback
        save_everystep=false, callback=callbacks, adaptive=false)

    return sol, semi

end

function set_inertia(time_method, source_terms, equations, problem_name)

    return ProblemSetup(problem_name=problem_name, equations=equations, initial_condition=initial_condition_gravity_wave,
        boundary_conditions=(x_neg=boundary_condition_periodic,
            x_pos=boundary_condition_periodic,
            y_neg=boundary_condition_slip_wall_2,
            y_pos=boundary_condition_slip_wall_2), tspan=(0.0, 1800.0), coordinates_min=(0.0, 0.0),
        coordinates_max=(300_000.0, 10_000.0), source_terms=source_terms, periodicity=(true, false), time_method=time_method)

end

function prim2perturb(u, x, equations::CompressibleEulerPotentialTemperatureEquations2DNC)
    rho, rho_v1, rho_v2, rho_theta, _ = u
    g = 9.81
    v2 = rho_v2 / rho
    v1 = rho_v1 / rho
    T0 = 250
    ps = 100000
    c_p = equations.c_p
    c_v = equations.c_v
    p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
    R = c_p - c_v
    delta = g / (R * T0)
    rhos = ps / (T0 * R)
    p0 = ps * exp(-delta * x[2])
    rho0 = rhos * exp(-delta * x[2])
    rhob = rho - rho0
    theta0 = (p0 / 100000)^(717 / 1004) * 100000 / (287 * rho0)
    pb = p - p0
    T = p / (rho * R)
    Tb = T - T0
    vb = v1 - 20.0
    thetab = (rho_theta - rhob * theta0 - rho0 * theta0) / rho0
    return vb, v2, pb, rhob, Tb, thetab
end

function prim2perturb(u, x, equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, rho_v1, rho_v2, rho_theta = u
    g = 9.81
    v2 = rho_v2 / rho
    v1 = rho_v1 / rho
    T0 = 250
    ps = 100000
    c_p = equations.c_p
    c_v = equations.c_v
    p = equations.p_0 * (equations.R * rho_theta / equations.p_0)^equations.gamma
    R = c_p - c_v
    delta = g / (R * T0)
    rhos = ps / (T0 * R)
    p0 = ps * exp(-delta * x[2])
    rho0 = rhos * exp(-delta * x[2])
    rhob = rho - rho0
    theta0 = (p0 / 100000)^(717 / 1004) * 100000 / (287 * rho0)
    pb = p - p0
    T = p / (rho * R)
    Tb = T - T0
    vb = v1 - 20.0
    thetab = (rho_theta - rhob * theta0 - rho0 * theta0) / rho0
    return vb, v2, pb, rhob, Tb, thetab
end

function prim2perturb(u, x, equations::CompressibleEulerEquations2DNC)
    rho, rho_v1, rho_v2, rho_e, _ = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    g = 9.81
    T0 = 250
    c_p = 1004.0
    c_v = 717.0
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
    R = c_p - c_v
    ps = 100_000
    T = p / (rho * R)
    Tb = T - T0
    delta = g / (R * T0)
    rhos = ps / (T0 * R)
    p0 = ps * exp(-delta * x[2])
    rho0 = rhos * exp(-delta * x[2])
    pb = p - p0
    rhob = rho - rho0
    rho_theta = (p / 100000)^(717 / 1004) * 100000 / 287
    theta0 = (p0 / 100000)^(717 / 1004) * 100000 / (287 * rho0)
    thetab = (rho_theta - rhob * theta0 - rho0 * theta0) / rho0
    vb = v1 - 20.0
    return vb, v2, pb, rhob, Tb, thetab
end

function prim2perturb(u, x, equations::CompressibleEulerEquations2D)
    rho, rho_v1, rho_v2, rho_e = u
    v1 = rho_v1 / rho
    v2 = rho_v2 / rho
    g = 9.81
    T0 = 250
    c_p = 1004.0
    c_v = 717.0
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
    R = c_p - c_v
    ps = 100_000
    T = p / (rho * R)
    Tb = T - T0
    delta = g / (R * T0)
    rhos = ps / (T0 * R)
    p0 = ps * exp(-delta * x[2])
    rho0 = rhos * exp(-delta * x[2])
    pb = p - p0
    rhob = rho - rho0
    rho_theta = (p / 100000)^(717 / 1004) * 100000 / 287
    theta0 = (p0 / 100000)^(717 / 1004) * 100000 / (287 * rho0)
    thetab = (rho_theta - rhob * theta0 - rho0 * theta0) / rho0
    vb = v1 - 20.0
    return vb, v2, pb, rhob, Tb, thetab
end


function get_perturbation!(sol, semi, equations, cells_per_dimension, polydeg)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u_wrap = Trixi.wrap_array(sol.u[end], semi)

    Tb = Array{Float64}(
        undef,
        polydeg + 1,
        polydeg + 1,
        cells_per_dimension[1] * cells_per_dimension[2],
    )

    wb = Array{Float64}(
        undef,
        polydeg + 1,
        polydeg + 1,
        cells_per_dimension[1] * cells_per_dimension[2],
    )

    ub = Array{Float64}(
        undef,
        polydeg + 1,
        polydeg + 1,
        cells_per_dimension[1] * cells_per_dimension[2],
    )

    pb = Array{Float64}(
        undef,
        polydeg + 1,
        polydeg + 1,
        cells_per_dimension[1] * cells_per_dimension[2],
    )

    thetab = Array{Float64}(
        undef,
        polydeg + 1,
        polydeg + 1,
        cells_per_dimension[1] * cells_per_dimension[2],
    )

    rhob = Array{Float64}(
        undef,
        polydeg + 1,
        polydeg + 1,
        cells_per_dimension[1] * cells_per_dimension[2],
    )

    @unpack node_coordinates = cache.elements
    for element in eachelement(solver, cache)
        for j in eachnode(solver), i in eachnode(solver)
            x_local =
                Trixi.get_node_coords(node_coordinates, equations, solver, i, j, element)
            ub[i, j, element], wb[i, j, element], pb[i, j, element], rhob[i, j, element], Tb[i, j, element], thetab[i, j, element] = prim2perturb(u_wrap[:, i, j, element], x_local, equations)
        end
    end
    return ub, wb, pb, rhob, Tb, thetab
end

function run_setups(; formulation, CFL, time_method, surface_flux, volume_flux, Nxmax)

    up_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    wp_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    rhop_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    Tp_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    pp_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    thetap_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])

    for polydeg in  (2, 3, 4)
        for n in  (1, 2, 4, 8)
            cells_per_dimension = (30 * n, 3 * n)
            @unpack equations = formulation
            sol, semi = inertia_gravity_wave(; polydeg=polydeg, CFL=CFL, time_method=time_method, cells_per_dimension=cells_per_dimension, surface_flux=surface_flux, volume_flux=volume_flux, problem_setup=formulation)
            ub, wb, pb, rhob, Tb, thetab = get_perturbation!(sol, semi, equations, cells_per_dimension, polydeg)
            ub_e, wb_e, pb_e, rhob_e, Tb_e, thetab_e = ComputeAnalyticalSolution(Nxmax, semi, cells_per_dimension, polydeg, 1800.0)
            l2, linf = compute_error_norms(wb, wb_e, semi.solver.basis.weights, semi)
            push!(wp_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(ub, ub_e, semi.solver.basis.weights, semi)
            push!(up_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(rhob, rhob_e, semi.solver.basis.weights, semi)
            push!(rhop_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(Tb, Tb_e, semi.solver.basis.weights, semi)
            push!(Tp_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(pb, pb_e, semi.solver.basis.weights, semi)
            push!(pp_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(thetab, thetab_e, semi.solver.basis.weights, semi)
            push!(thetap_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
        end
    end
    if isa(volume_flux, Tuple)
        stringflux = string(volume_flux[1])
        stringfluxnc = string(volume_flux[2])
    else
        stringflux = string(volume_flux)
        stringfluxnc = "_"
    end
    CSV.write(pwd() * "/test_cases/gravitywaves/out/up_error_CFL_"*string(CFL) * formulation.problem_name * "_" * stringflux * "_" * stringfluxnc *"_.csv", up_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/wp_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", wp_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/pp_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", pp_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/Tp_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", Tp_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/rhop_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", rhop_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/thetap_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", thetap_error)
    return nothing
end

function run_setups_warped(; formulation, CFL, time_method, surface_flux, volume_flux, Nxmax)

    up_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    wp_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    rhop_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    Tp_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    pp_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])
    thetap_error = DataFrame(polydeg=Int[], problem_name=String[], dx=Float64[], l2_error=Float64[], linf_error=Float64[])

    for polydeg in  (2, 3, 4)
        for n in  (1, 2, 4, 8, 10)
            cells_per_dimension = (30 * n, 3 * n)
            @unpack equations = formulation
            sol, semi = inertia_gravity_wave_warped(; polydeg=polydeg, CFL=CFL, time_method=time_method, cells_per_dimension=cells_per_dimension, surface_flux=surface_flux, volume_flux=volume_flux, problem_setup=formulation)
            # pd = PlotData2D(sol)
            # a = Plots.plot(getmesh(pd), aspect_ratio = 10)
            # Plots.savefig(a, "mesh.pdf")
            # throw(error)
            ub, wb, pb, rhob, Tb, thetab = get_perturbation!(sol, semi, equations, cells_per_dimension, polydeg)
            ub_e, wb_e, pb_e, rhob_e, Tb_e, thetab_e = ComputeAnalyticalSolution(Nxmax, semi, cells_per_dimension, polydeg, 1800.0)
            l2, linf = compute_error_norms(wb, wb_e, semi.solver.basis.weights, semi)
            push!(wp_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(ub, ub_e, semi.solver.basis.weights, semi)
            push!(up_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(rhob, rhob_e, semi.solver.basis.weights, semi)
            push!(rhop_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(Tb, Tb_e, semi.solver.basis.weights, semi)
            push!(Tp_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(pb, pb_e, semi.solver.basis.weights, semi)
            push!(pp_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
            l2, linf = compute_error_norms(thetab, thetab_e, semi.solver.basis.weights, semi)
            push!(thetap_error, (polydeg, formulation.problem_name, 10000.0 / n, l2, linf))
        end
    end
    if isa(volume_flux, Tuple)
        stringflux = string(volume_flux[1])
        stringfluxnc = string(volume_flux[2])
    else
        stringflux = string(volume_flux)
        stringfluxnc = "_"
    end
    CSV.write(pwd() * "/test_cases/gravitywaves/out/warped_up_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux* "_" * stringfluxnc * "_.csv", up_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/warped_wp_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", wp_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/warped_pp_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", pp_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/warped_Tp_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux * "_" * stringfluxnc *"_.csv", Tp_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/warped_rhop_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", rhop_error)
    CSV.write(pwd() * "/test_cases/gravitywaves/out/warped_thetap_error_CFL_"*string(CFL)  * formulation.problem_name * "_" * stringflux *"_" * stringfluxnc * "_.csv", thetap_error)
    return nothing
end


function run_analytical_solution(; CFL = 0.1, flux_nonconservative = flux_nonconservative_gravity_log, flux_pot = flux_theta)

    time_method = SSPRK43(thread = Trixi.True())
    Nxmax = 1200
    cs = sqrt(1004.0 / 717.0 * 287.0 * 250.0)

    PT = set_inertia(time_method, source_terms_gravity, CompressibleEulerPotentialTemperatureEquations2D(), "PT")

    surface_flux = FluxLMARS(cs)
    volume_flux = flux_pot
    run_setups(; formulation=PT, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax=Nxmax)

    PT_NC = set_inertia(time_method, nothing, CompressibleEulerPotentialTemperatureEquations2DNC(), "PT_NC")

    surface_flux = (FluxLMARS(cs), flux_nonconservative)
    volume_flux = (flux_pot, flux_nonconservative)

    run_setups(; formulation=PT_NC, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax)

    EU_NC = set_inertia(time_method, nothing, CompressibleEulerEquations2DNC(1004.0 / 717.0), "EUNC")

    surface_flux = (FluxLMARS(cs), flux_nonconservative)
    volume_flux = (flux_ranocha, flux_nonconservative)

    run_setups(; formulation=EU_NC, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax)

    EU = set_inertia(time_method, source_terms_gravity, CompressibleEulerEquations2D(1004.0 / 717.0), "Euler")

    surface_flux = FluxLMARS(cs)
    volume_flux = flux_ranocha
    run_setups(; formulation=EU, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax)

    return nothing
end
@inline function flux_zero(u_ll, u_rr, normal_direction, equations)
    return zero(u_ll)
end

function run_analytical_solution_warped(; CFL = 0.1, flux_nonconservative = flux_nonconservative_gravity_log)

    time_method = SSPRK43(thread = Trixi.True())
    Nxmax = 1200
    cs = sqrt(1004.0 / 717.0 * 287.0 * 250.0)

    PT = set_inertia(time_method, source_terms_gravity, CompressibleEulerPotentialTemperatureEquations2D(), "PT")

    surface_flux = FluxLMARS(cs)
    volume_flux = flux_theta
    run_setups_warped(; formulation=PT, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax=Nxmax)

    PT_NC = set_inertia(time_method, nothing, CompressibleEulerPotentialTemperatureEquations2DNC(), "PT_NC_zero")

    surface_flux = (FluxLMARS(cs), flux_zero)
    volume_flux = (flux_theta, flux_nonconservative)

    run_setups_warped(; formulation=PT_NC, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax)

    EU_NC = set_inertia(time_method, nothing, CompressibleEulerEquations2DNC(1004.0 / 717.0), "EUNC_zero")

    surface_flux = (FluxLMARS(cs), flux_zero)
    volume_flux = (flux_ranocha, flux_nonconservative)

    run_setups_warped(; formulation=EU_NC, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax)

    EU = set_inertia(time_method, source_terms_gravity, CompressibleEulerEquations2D(1004.0 / 717.0), "Euler")

    surface_flux = FluxLMARS(cs)
    volume_flux = flux_ranocha
    run_setups_warped(; formulation=EU, CFL=CFL, time_method=time_method, surface_flux=surface_flux, volume_flux=volume_flux, Nxmax)

    return nothing
end

run_analytical_solution()
run_analytical_solution_warped()
