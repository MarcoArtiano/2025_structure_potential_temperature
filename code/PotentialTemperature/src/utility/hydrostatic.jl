function prim2velocity(u, x, equations::Union{CompressibleEulerEquations1D,CompressibleEulerPotentialTemperatureEquations1D}) 
    rho, rho_v1, _ = u

    v1 = rho_v1 / rho

    return v1
end

function prim2velocity(u, x, equations::Union{CompressibleEulerEquations1DNC,CompressibleEulerPotentialTemperatureEquations1DNC}) 
    rho, rho_v1, _, _ = u

    v1 = rho_v1 / rho

    return v1
end

function get_velocity(sol, semi, equations::Union{CompressibleEulerEquations1D,CompressibleEulerEquations1DNC,CompressibleEulerPotentialTemperatureEquations1D, CompressibleEulerPotentialTemperatureEquations1DNC}, initial_refinement_level, polydeg)
    mesh, equations, solver, cache = Trixi.mesh_equations_solver_cache(semi)
    u_wrap = Trixi.wrap_array(sol.u, semi)
    velocity = Array{Float64}(
        undef,
        polydeg + 1,
        2^initial_refinement_level,
    )
    @unpack node_coordinates = cache.elements
    for element in eachelement(solver, cache)
        for i in eachnode(solver)
            x_local =
                Trixi.get_node_coords(node_coordinates, equations, solver, i, element)
            velocity[i, element] =
                prim2velocity(u_wrap[:, i, element], x_local, equations)
        end
    end

    return velocity

end