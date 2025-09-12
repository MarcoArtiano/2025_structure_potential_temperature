module PotentialTemperature

using MuladdMacro

# Conservative Formulation
include("equations/compressible_euler_potential_temperature_1d.jl")
include("equations/compressible_euler_potential_temperature_2d.jl")
include("equations/compressible_euler_potential_temperature_3d.jl")

# Gravity-non Conservative Formulation
include("equations/compressible_euler_1d_nc.jl")
include("equations/compressible_euler_2d_nc.jl")
include("equations/compressible_euler_3d_nc.jl")
include("equations/compressible_euler_potential_temperature_1d_nc.jl")
include("equations/compressible_euler_potential_temperature_2d_nc.jl")
include("equations/compressible_euler_potential_temperature_3d_nc.jl")

export CompressibleEulerEquations1DNC,
    CompressibleEulerEquations2DNC,
    CompressibleEulerEquations3DNC,
    CompressibleEulerPotentialTemperatureEquations1D,
    CompressibleEulerPotentialTemperatureEquations2D,
    CompressibleEulerPotentialTemperatureEquations3D,
    CompressibleEulerPotentialTemperatureEquations1DNC,
    CompressibleEulerPotentialTemperatureEquations2DNC,
    CompressibleEulerPotentialTemperatureEquations3DNC

export flux_theta, flux_theta_AM, flux_theta_rhos, flux_theta_rhos_AM, flux_theta_global, flux_theta_sto, flux_theta_es
export flux_theta_mod
export pressure, total_energy, energy_kinetic, entropy

export flux_nonconservative_gravity_log, flux_nonconservative_gravity_gamma, flux_nonconservative_gravity_am
export flux_nonconservative_gravity_gamma_curvilinear
export well_balanced_v1, well_balanced_v2, flux_nonconservative_gravity_gamma_wb, flux_nonconservative_gravity_gamma_surface
export flux_theta_surface
## Optimized Stolarsky Mean for Potential Temperature, gaining almost a factor of two.
@inline function stolarsky_mean_opt(x::RealT, y::RealT, gamma::RealT, xg::RealT, yg::RealT, equations) where {RealT<:Real}
    epsilon_f2 = convert(RealT, 1.0e-4)
    f2 = (x * (x - 2 * y) + y * y) / (x * (x + 2 * y) + y * y) # f2 = f^2
    if f2 < epsilon_f2
        # convenience coefficients
        c1 = convert(RealT, 1 / 3) * (gamma - 2)
        c2 = convert(RealT, -1 / 15) * (gamma + 1) * (gamma - 3) * c1
        c3 = convert(RealT, -1 / 21) * (2 * gamma * (gamma - 2) - 9) * c2
        return 0.5f0 * (x + y) * @evalpoly(f2, 1, c1, c2, c3)
    else
        return equations.stolarsky_factor * (yg - xg) * x * y /
               (yg * x - xg * y)
    end
end

@inline function source_terms_gravity(u, x, t,
    equations::CompressibleEulerEquations2D)
    rho, _, rho_v2, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -9.81 * rho,
        -9.81 * rho_v2)
end

@inline function source_terms_gravity(u, x, t,
    equations::Union{CompressibleEulerEquations1D,CompressibleEulerPotentialTemperatureEquations1D})
    rho, rho_v1, _ = u
    return SVector(zero(eltype(u)), -9.81 * rho,
        -9.81 * rho_v1)
end

@inline function source_terms_gravity(u, x, t,
    equations::CompressibleEulerPotentialTemperatureEquations2D)
    rho, _, _, _ = u
    return SVector(zero(eltype(u)), zero(eltype(u)), -equations.g * rho,
        zero(eltype(u)))
end

Base.@kwdef struct ProblemSetup{NameT,EquationsT,InitialConditionT,BoundaryConditionT,TspanT,CoordinateT,SourceT,PeriodicityT,TimeMethodT}
    problem_name::NameT
    equations::EquationsT
    initial_condition::InitialConditionT
    boundary_conditions::BoundaryConditionT
    tspan::TspanT
    coordinates_min::CoordinateT
    coordinates_max::CoordinateT
    source_terms::SourceT
    periodicity::PeriodicityT
    time_method::TimeMethodT
end

export ProblemSetup


@inline function boundary_condition_slip_wall_2(u_inner, orientation_or_normal, direction,
    x, t,
    surface_flux_function,
    equations::Union{CompressibleEulerEquations1D,CompressibleEulerPotentialTemperatureEquations1D})

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
        -u_inner[2],
        u_inner[3])

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

@inline function boundary_condition_slip_wall_2(u_inner, normal_direction::AbstractVector,
    x, t,
    surface_flux_function,
    equations::Union{CompressibleEulerEquations3D,CompressibleEulerPotentialTemperatureEquations3D})
    # normalize the outward pointing direction
    normal = normal_direction / Trixi.norm(normal_direction)
    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3] + normal[3] * u_inner[4]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
        u_inner[2] - 2 * u_normal * normal[1],
        u_inner[3] - 2 * u_normal * normal[2],
        u_inner[4] - 2 * u_normal * normal[3],
        u_inner[5])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

    return flux
end

@inline function boundary_condition_slip_wall_2(u_inner, orientation, direction, x, t,
    surface_flux_function,
    equations::Union{CompressibleEulerEquations3D,CompressibleEulerPotentialTemperatureEquations3D})
    # Boundary state is equal to the inner state except for the velocity. For boundaries
    # in the -x/+x direction, we multiply the velocity in the x direction by -1.
    # Similarly, for boundaries in the -y/+y or -z/+z direction, we multiply the
    # velocity in the y or z direction by -1
    if direction in (1, 2) # x direction
        u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4],
            u_inner[5])
    elseif direction in (3, 4) # y direction
        u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4],
            u_inner[5])
    else # z direction = (5, 6)
        u_boundary = SVector(u_inner[1], u_inner[2], u_inner[3], -u_inner[4],
            u_inner[5])
    end

    # Calculate boundary flux depending on the orientation of the boundary
    # Odd directions are in negative coordinate direction, 
    # even directions are in positive coordinate direction.
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end

@inline function boundary_condition_slip_wall_2(u_inner, normal_direction::AbstractVector,
    x, t,
    surface_flux_function,
    equations::Union{CompressibleEulerPotentialTemperatureEquations2D, CompressibleEulerEquations2D})
    # normalize the outward pointing direction
    normal = normal_direction / norm(normal_direction)

    # compute the normal velocity
    u_normal = normal[1] * u_inner[2] + normal[2] * u_inner[3]

    # create the "external" boundary solution state
    u_boundary = SVector(u_inner[1],
        u_inner[2] - 2 * u_normal * normal[1],
        u_inner[3] - 2 * u_normal * normal[2],
        u_inner[4])

    # calculate the boundary flux
    flux = surface_flux_function(u_inner, u_boundary, normal_direction, equations)

    return flux
end

"""
boundary_condition_slip_wall(u_inner, orientation, direction, x, t,
surface_flux_function, equations::ShallowWaterEquations2D)

Should be used together with [`TreeMesh`](@ref).
"""
@inline function boundary_condition_slip_wall_2(u_inner, orientation,
    direction, x, t,
    surface_flux_function,
    equations::Union{CompressibleEulerPotentialTemperatureEquations2D, CompressibleEulerEquations2D})
    ## get the appropriate normal vector from the orientation
    if orientation == 1
        u_boundary = SVector(u_inner[1], -u_inner[2], u_inner[3], u_inner[4])
    else # orientation == 2
        u_boundary = SVector(u_inner[1], u_inner[2], -u_inner[3], u_inner[4])
    end

    # Calculate boundary flux
    if iseven(direction) # u_inner is "left" of boundary, u_boundary is "right" of boundary
        flux = surface_flux_function(u_inner, u_boundary, orientation, equations)
    else # u_boundary is "left" of boundary, u_inner is "right" of boundary
        flux = surface_flux_function(u_boundary, u_inner, orientation, equations)
    end

    return flux
end

export source_terms_gravity, boundary_condition_slip_wall_2

include("utility/plots.jl")
export load_processed_data, create_plot!, create_plot2!, plot_var!

include("utility/hydrostatic.jl")
export prim2velocity, get_velocity

include("utility/baroclinic.jl")
export perturbation_stream_function, basic_state_baroclinic_instability_longitudinal_velocity, cartesian_to_sphere

include("utility/contours.jl")
export ContourData

end # module PotentialTemperature
