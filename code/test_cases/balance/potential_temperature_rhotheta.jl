using OrdinaryDiffEqSSPRK
using PotentialTemperature.Trixi
using PotentialTemperature
using CairoMakie
using DoubleFloats
using LaTeXStrings
using DelimitedFiles

function initial_condition_adiabatic(x, t, equations::Union{CompressibleEulerPotentialTemperatureEquations2DNC,CompressibleEulerEquations2DNC})
    RealT = eltype(x)
    g = RealT(9.81)
    c_p = RealT(1004.0)
    c_v = RealT(717.0)
    # center of perturbation
    p0 = RealT(100_000)
    # perturbation in potential temperature
    R = c_p - c_v    # gas constant (dry air)
    potential_temperature = 300.0

    # Exner pressure, solves hydrostatic equation for x[2]
    exner = RealT(1) - g / (c_p * potential_temperature) * x[2]

    # pressure
    R = c_p - c_v    # gas constant (dry air)
    p = p0 * exner^(c_p / R)

    # temperature
    T = potential_temperature * exner
    # density
    rho = p / (R * T)
    v1 = RealT(0.0)
    v2 = RealT(0.0)

    return prim2cons(SVector(rho, v1, v2, p, x[2]), equations)
end

@inline function flux_zero(u_ll, u_rr, normal_direction, equations)
    return zero(u_ll)
end

RealT = Double64
T = RealT(5000)
equations = CompressibleEulerPotentialTemperatureEquations2DNC()
polydeg = 2
basis = LobattoLegendreBasis(RealT, polydeg)
cs = RealT(340.0)
flux_nonconservative = flux_nonconservative_gravity_gamma
surface_flux = (FluxLMARS(cs), flux_zero)
volume_flux = (flux_theta, flux_nonconservative)
volume_integral = VolumeIntegralFluxDifferencing(volume_flux)
solver = DGSEM(basis, surface_flux, volume_integral)

trees_per_dimension = (16, 16)

function mapping(xi, eta)
    x = xi + RealT(0.1) * sinpi(xi) * sinpi(eta)
    y = eta + RealT(0.1) * sinpi(xi) * sinpi(eta)
    return SVector(RealT(1000) * RealT(0.5) * (RealT(1) + x), RealT(1000) * RealT(0.5) * (RealT(1) + y))
end

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg, 
mapping=mapping, periodicity = (false,false), initial_refinement_level=0, RealT = RealT)

boundary_conditions =  Dict( :x_pos => boundary_condition_slip_wall_2, 
                             :x_neg => boundary_condition_slip_wall_2,
                             :y_pos => boundary_condition_slip_wall_2,
                             :y_neg => boundary_condition_slip_wall_2)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_adiabatic, solver, boundary_conditions=boundary_conditions)

dt = RealT(0.01)
tspan = (zero(RealT), T)

summary_callback = SummaryCallback()

analysis_interval = 1
analysis_callback = AnalysisCallback(semi, interval=analysis_interval, extra_analysis_integrals=(well_balanced_v1, well_balanced_v2), save_analysis = true, output_directory = pwd()*"/test_cases/balance/out", analysis_filename = "well_balancing_theta_$(RealT)_rhotheta.dat")

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

ode = semidiscretize(semi, tspan)

sol = solve(ode, SSPRK43(); dt = dt, ode_default_options()..., callback = callbacks, adaptive = false)

summary_callback();

RealT = Float64
equations = CompressibleEulerPotentialTemperatureEquations2DNC()
basis = LobattoLegendreBasis(RealT, polydeg)
cs = RealT(340.0)
surface_flux = (FluxLMARS(cs), flux_nonconservative)

solver = DGSEM(basis, surface_flux, volume_integral)

mesh = P4estMesh(trees_per_dimension, polydeg = polydeg, mapping=mapping, periodicity = (false,false), initial_refinement_level=0, RealT = RealT)

boundary_conditions =  Dict( :x_pos => boundary_condition_slip_wall_2, 
                             :x_neg => boundary_condition_slip_wall_2,
                             :y_pos => boundary_condition_slip_wall_2,
                             :y_neg => boundary_condition_slip_wall_2)

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_adiabatic, solver, boundary_conditions=boundary_conditions)

dt = RealT(0.01)
tspan = (zero(RealT), T)
ode = semidiscretize(semi, tspan)

summary_callback = SummaryCallback()

analysis_interval = 1
analysis_callback = AnalysisCallback(semi, interval=analysis_interval, extra_analysis_integrals=(well_balanced_v1, well_balanced_v2), save_analysis = true, output_directory = pwd()*"/test_cases/balance/out", analysis_filename = "well_balancing_theta_$(RealT)_rhotheta.dat")

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

sol = solve(ode, SSPRK43(); dt = dt, ode_default_options()..., callback = callbacks, adaptive = false)

summary_callback();
