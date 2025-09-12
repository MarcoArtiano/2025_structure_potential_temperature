using OrdinaryDiffEqSSPRK
using PotentialTemperature.Trixi
using PotentialTemperature
using CairoMakie
using DoubleFloats
using LaTeXStrings
using DelimitedFiles

function initial_condition_isothermal(x, t, equations::Union{CompressibleEulerPotentialTemperatureEquations2DNC,CompressibleEulerEquations2DNC})
    RealT = eltype(x)
    g = RealT(9.81)
    c_p = RealT(1004.0)
    c_v = RealT(717.0)
    # center of perturbation
    T0 = RealT(250.0)
    p0 = RealT(100_000)
    # perturbation in potential temperature
    R = c_p - c_v    # gas constant (dry air)
   delta = g / (R * T0)

    rho0 = p0 / (T0 * R)
    p = p0 * exp(-delta * x[2])
    rho = rho0 * exp(-delta * x[2])
    v1 = zero(RealT)
    v2 = zero(RealT)

    #return prim2cons(SVector(rho, v1, v2, p, g * x[2]), equations)
    return prim2cons(SVector(rho, v1, v2, p, x[2]), equations)

end

@inline function flux_zero(u_ll, u_rr, normal_direction, equations)
    return zero(u_ll)
end

RealT = Float64
T = RealT(0.0)
equations = CompressibleEulerEquations2DNC(RealT(1004/717))
polydeg = 2
basis = LobattoLegendreBasis(RealT, polydeg)
cs = RealT(340.0)
flux_nonconservative = flux_nonconservative_gravity_log
surface_flux = (FluxLMARS(cs), flux_nonconservative)
volume_flux = (flux_ranocha, flux_nonconservative)
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

semi = SemidiscretizationHyperbolic(mesh, equations, initial_condition_isothermal, solver, boundary_conditions=boundary_conditions)

dt = RealT(0.01)
tspan = (zero(RealT), T)

summary_callback = SummaryCallback()

analysis_interval = 1
analysis_callback = AnalysisCallback(semi, interval=analysis_interval, extra_analysis_integrals=(well_balanced_v1, well_balanced_v2), save_analysis = true, output_directory = pwd()*"/test_cases/balance/out", analysis_filename = "well_balancing_$(RealT)_mesh.dat")

alive_callback = AliveCallback(analysis_interval = analysis_interval)

callbacks = CallbackSet(summary_callback, analysis_callback, alive_callback)

ode = semidiscretize(semi, tspan)

sol = solve(ode, SSPRK43(); dt = dt, ode_default_options()..., callback = callbacks, adaptive = false)

summary_callback();

x, y, data = ContourData(sol.u[end], semi, trees_per_dimension, equations)

h = Figure(size = (850, 400))
labelsize = 20
kwargs = (xlabel = L"$x$ [km]", xlabelsize = labelsize, ylabelsize = labelsize, limits = ((0, 1), (0, 1)), xticklabelsize = 17.0, yticklabelsize = 17.0)

Axis(h[1, 1]; kwargs..., ylabel = L"$z$ [km]", aspect = DataAspect())

for j in axes(y, 2)
    lines!(x[:, j] ./ 1e3, y[:, j] ./ 1e3, color = :gray, linewidth = 0.5)
end

for i in axes(x, 1)
    lines!(x[i, :] ./ 1e3, y[i, :] ./ 1e3, color = :gray, linewidth = 0.5)
end

save(pwd() * "/test_cases/balance/plots/mesh.pdf", h)