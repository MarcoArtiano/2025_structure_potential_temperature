using OrdinaryDiffEqSSPRK
using PotentialTemperature.Trixi
using LinearAlgebra
using PotentialTemperature
###############################################################################
# Setup for the baroclinic instability test
include("plots.jl")

function initial_condition_baroclinic_instability(x, t,
    equations::Union{CompressibleEulerPotentialTemperatureEquations3D,CompressibleEulerEquations3D})
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = 6.371229e6
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Stream function type perturbation
    u_perturbation, v_perturbation = perturbation_stream_function(lon, lat, z)

    u += u_perturbation
    v = v_perturbation

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u - sin(lat) * cos(lon) * v
    v2 = cos(lon) * u - sin(lat) * sin(lon) * v
    v3 = cos(lat) * v

    return prim2cons(SVector(rho, v1, v2, v3, p), equations)
end
# Initial condition for an idealized baroclinic instability test
# https://doi.org/10.1002/qj.2241, Section 3.2 and Appendix A
function initial_condition_baroclinic_instability(x, t,
    equations::Union{CompressibleEulerPotentialTemperatureEquations3DNC,CompressibleEulerEquations3DNC})
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = 6.371229e6
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Stream function type perturbation
    u_perturbation, v_perturbation = perturbation_stream_function(lon, lat, z)

    u += u_perturbation
    v = v_perturbation

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u - sin(lat) * cos(lon) * v
    v2 = cos(lon) * u - sin(lat) * sin(lon) * v
    v3 = cos(lat) * v
    radius_earth = 6.371229e6  # a
    gravitational_acceleration = 9.81    # g

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)
    if z > 0
        r = norm(x)
    else
        r = -(2 * radius_earth^3) / (x[1]^2 + x[2]^2 + x[3]^2)
    end
    r = -norm(x)
    phi = radius_earth^2 * gravitational_acceleration / r

    return prim2cons(SVector(rho, v1, v2, v3, p, phi), equations)
end

# Steady state for RHS correction below
function steady_state_baroclinic_instability(x, t, equations::Union{CompressibleEulerEquations3D,CompressibleEulerPotentialTemperatureEquations3D})
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = 6.371229e6
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u
    v2 = cos(lon) * u
    v3 = 0.0

    return prim2cons(SVector(rho, v1, v2, v3, p), equations)
end
# Steady state for RHS correction below
function steady_state_baroclinic_instability(x, t, equations::Union{CompressibleEulerPotentialTemperatureEquations3DNC,CompressibleEulerEquations3DNC})
    lon, lat, r = cartesian_to_sphere(x)
    radius_earth = 6.371229e6
    # Make sure that the r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)

    # Unperturbed basic state
    rho, u, p = basic_state_baroclinic_instability_longitudinal_velocity(lon, lat, z)

    # Convert spherical velocity to Cartesian
    v1 = -sin(lon) * u
    v2 = cos(lon) * u
    v3 = 0.0
    radius_earth = 6.371229e6  # a
    gravitational_acceleration = 9.81     # g

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)

    if z > 0
        r = norm(x)
    else
        r = -(2 * radius_earth^3) / (x[1]^2 + x[2]^2 + x[3]^2)
    end
    r = -norm(x)
    phi = radius_earth^2 * gravitational_acceleration / r

    return prim2cons(SVector(rho, v1, v2, v3, p, phi), equations)
end



@inline function source_terms_baroclinic_instability(u, x, t,
    equations::Union{CompressibleEulerPotentialTemperatureEquations3DNC,CompressibleEulerEquations3DNC})
    radius_earth = 6.371229e6  # a
    gravitational_acceleration = 9.81     # g
    angular_velocity = 7.29212e-5  # Ω

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)
    r = z + radius_earth

    du1 = zero(eltype(u))
    du2 = zero(eltype(u))
    du3 = zero(eltype(u))
    du4 = zero(eltype(u))
    du5 = zero(eltype(u))
    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 -= -2 * angular_velocity * u[3]
    du3 -= 2 * angular_velocity * u[2]

    return SVector(du1, du2, du3, du4, du5, zero(eltype(u)))
end

@inline function source_terms_baroclinic_instability(u, x, t,
    equations::CompressibleEulerEquations3D)
    radius_earth = 6.371229e6  # a
    gravitational_acceleration = 9.81     # g
    angular_velocity = 7.29212e-5  # Ω

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)
    r = z + radius_earth

    du1 = zero(eltype(u))

    # Gravity term
    temp = -gravitational_acceleration * radius_earth^2 / r^3
    du2 = temp * u[1] * x[1]
    du3 = temp * u[1] * x[2]
    du4 = temp * u[1] * x[3]
    du5 = temp * (u[2] * x[1] + u[3] * x[2] + u[4] * x[3])

    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 -= -2 * angular_velocity * u[3]
    du3 -= 2 * angular_velocity * u[2]

    return SVector(du1, du2, du3, du4, du5)
end

@inline function source_terms_baroclinic_instability(u, x, t,
    equations::CompressibleEulerPotentialTemperatureEquations3D)
    radius_earth = 6.371229e6  # a
    gravitational_acceleration = 9.81     # g
    angular_velocity = 7.29212e-5  # Ω

    r = norm(x)
    # Make sure that r is not smaller than radius_earth
    z = max(r - radius_earth, 0.0)
    r = z + radius_earth

    du1 = zero(eltype(u))

    # Gravity term
    temp = -gravitational_acceleration * radius_earth^2 / r^3
    du2 = temp * u[1] * x[1]
    du3 = temp * u[1] * x[2]
    du4 = temp * u[1] * x[3]
    du5 = zero(eltype(u))

    # Coriolis term, -2Ω × ρv = -2 * angular_velocity * (0, 0, 1) × u[2:4]
    du2 -= -2 * angular_velocity * u[3]
    du3 -= 2 * angular_velocity * u[2]

    return SVector(du1, du2, du3, du4, du5)
end

###############################################################################
# Start of the actual elixir, semidiscretization of the problem

function main(equations, surface_flux, volume_flux, T, filename, trees_per_cube_face, polydeg)

    initial_condition = initial_condition_baroclinic_instability

    boundary_conditions = Dict(:inside => boundary_condition_slip_wall_2,
        :outside => boundary_condition_slip_wall_2)
        
    solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux,
        volume_integral=VolumeIntegralFluxDifferencing(volume_flux))

    mesh = Trixi.P4estMeshCubedSphere(trees_per_cube_face..., 6.371229e6, 30000.0,
        polydeg=polydeg, initial_refinement_level=0)

    semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition, solver,
        source_terms=source_terms_baroclinic_instability,
        boundary_conditions=boundary_conditions)

    ###############################################################################
    # ODE solvers, callbacks etc.

    tspan = (0.0, T * 24 * 60 * 60.0) # time in seconds for 10 days

    ode = semidiscretize(semi, tspan)

    summary_callback = SummaryCallback()

    analysis_interval = 10000
    analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

    alive_callback = AliveCallback(analysis_interval=analysis_interval)

    save_solution = SaveSolutionCallback(interval=300000, save_initial_solution=true,
        save_final_solution=true,
        output_directory="test_cases/baroclinic/" * filename * "$(trees_per_cube_face[1])x$(trees_per_cube_face[2])_$(polydeg)")

    callbacks = CallbackSet(summary_callback,
        analysis_callback,
        alive_callback,
        save_solution)

    tol = 1e-6
    ###############################################################################
    # Use a Runge-Kutta method with automatic (error based) time step size control
    # Enable threading of the RK method for better performance on multiple threads

    sol = solve(ode,
        SSPRK43(thread = Trixi.True());
        abstol=tol,
        reltol=tol, ode_default_options()...,
        callback=callbacks)
        return sol,semi
end

function run_baroclinic(; T, trees_per_cube_face=(8, 4), polydeg = 5, flux_nonconservative = flux_nonconservative_gravity_log)
     folder = pwd() * "/test_cases/baroclinic/data/"
     
     volume_flux = flux_ranocha
     surface_flux = FluxLMARS(340.0)
     equations = CompressibleEulerEquations3D(1004 / 717)
     filename = "Euler_" * string(T)
     sol, semi = main(equations, surface_flux, volume_flux, T, filename, trees_per_cube_face, polydeg)

     surface_flux = (FluxLMARS(340), flux_nonconservative)
     volume_flux = (flux_ranocha, flux_nonconservative)
     equations = CompressibleEulerEquations3DNC(1004 / 717)
     equations_euler = equations
     filename = "NCEuler_" * string(T)
     sol_euler, semi_euler =  main(equations, surface_flux, volume_flux, T, filename, trees_per_cube_face, polydeg)

     volume_flux = (flux_theta, flux_nonconservative)
     equations = CompressibleEulerPotentialTemperatureEquations3DNC()
     equations_theta = equations
     filename = "NCPotential_" * string(T)
     sol_theta, semi_theta = main(equations, surface_flux, volume_flux, T, filename, trees_per_cube_face, polydeg)
     include("test_cases/baroclinic/plots.jl")
     contour_baroclinic(sol_euler, semi_euler,sol_theta, semi_theta, trees_per_cube_face, 2* (polydeg+  1), equations_euler, equations_theta)
end

function run_single(; T, trees_per_cube_face=(8,4), polydeg = 5)
    volume_flux = flux_kennedy_gruber
    surface_flux = FluxLMARS(340.0)
    equations = CompressibleEulerEquations3D(1004 / 717)
    filename = "Euler" * string(T)
    sol, semi = main(equations, surface_flux, volume_flux, T, filename, trees_per_cube_face, polydeg)
    return sol,semi
end

run_baroclinic(T = 10.0)