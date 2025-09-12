function ContourData(sol, semi, cells_per_dimension, equations)
    
    node_coordinates = semi.cache.elements.node_coordinates
    polydeg = Trixi.nnodes(semi.solver.basis) - 1
    Nx = polydeg * cells_per_dimension[1] + 1
    Ny = polydeg * cells_per_dimension[2] + 1
    u = Trixi.wrap_array(sol, semi)

    x = zeros(Ny, Nx)
    y = copy(x)
    # u = sol.u
    nvars = size(u,1)
    @show Ny, Nx
    data = zeros(Float64, nvars+2, Ny, Nx)
    
    linear_indices = LinearIndices(cells_per_dimension)
    for cell_x in 1:cells_per_dimension[1]
        for cell_y in 1:cells_per_dimension[2]
            element = linear_indices[cell_x, cell_y]
            
            for jlocal in 1:polydeg
                for ilocal in 1:polydeg
                    i, j = compute_global_index(element, ilocal, jlocal, polydeg, cell_x, cell_y)
                    x[j, i] = node_coordinates[1, ilocal, jlocal, element]
                    y[j, i] = node_coordinates[2, ilocal, jlocal, element]
                    data[1:nvars, j, i] = u[:, ilocal, jlocal, element]
                    data[end-1,j,i] = cons2theta(u[:,ilocal,jlocal,element], node_coordinates[2, ilocal, jlocal, element], equations)
                    data[end, j, i] = element
                end
            end

            if cell_x == cells_per_dimension[1]
                ilocal = polydeg + 1
                for jlocal in 1:polydeg
                    i, j = compute_global_index(element, ilocal, jlocal, polydeg, cell_x, cell_y)
                    x[j, i] = node_coordinates[1, ilocal, jlocal, element]
                    y[j, i] = node_coordinates[2, ilocal, jlocal, element]
                    data[1:nvars, j, i] = u[:, ilocal, jlocal, element]
                    data[end-1,j,i] = cons2theta(u[:,ilocal,jlocal,element], node_coordinates[2, ilocal, jlocal, element], equations)
                    data[end, j, i] = element
                end
            end

            if cell_y == cells_per_dimension[2]
                jlocal = polydeg + 1
                for ilocal in 1:polydeg
                    i, j = compute_global_index(element, ilocal, jlocal, polydeg, cell_x, cell_y)
                    x[j, i] = node_coordinates[1, ilocal, jlocal, element]
                    y[j, i] = node_coordinates[2, ilocal, jlocal, element]
                    data[1:nvars, j, i] = u[:, ilocal, jlocal, element]
                    data[end-1,j,i] = cons2theta(u[:,ilocal,jlocal,element],node_coordinates[2, ilocal, jlocal, element],equations)
                    data[end, j, i] = element
                end
            end

            if cell_y == cells_per_dimension[2] && cell_x == cells_per_dimension[1]
                ilocal = polydeg + 1
                jlocal = polydeg + 1
                    i, j = compute_global_index(element, ilocal, jlocal, polydeg, cell_x, cell_y)
                    x[j, i] = node_coordinates[1, ilocal, jlocal, element]
                    y[j, i] = node_coordinates[2, ilocal, jlocal, element]
                    data[1:nvars, j, i] = u[:, ilocal, jlocal, element]
                    data[end-1,j,i] = cons2theta(u[:,ilocal,jlocal,element], node_coordinates[2, ilocal, jlocal, element], equations)
                    data[end, j, i] = element
            end
        end
    end
    return x, y, data
end

function cons2theta(u,z, equations::Union{CompressibleEulerPotentialTemperatureEquations2D,CompressibleEulerPotentialTemperatureEquations2DNC})
    rho, rho_v1, rho_v2, rho_theta = u

    theta = rho_theta/rho - 300 * exp(0.01^2 / equations.g * z)

    return theta
end

function cons2theta(u, z, equations::Union{CompressibleEulerEquations2D,CompressibleEulerEquations2DNC})
	rho, rho_v1, rho_v2, rho_e = u
    v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	c_p = 1004.0
	c_v = 717.0
    p = (equations.gamma - 1) * (rho_e - 0.5 * (rho_v1 * v1 + rho_v2 * v2))
	R = c_p - c_v
	p0 = 100000
	rho_theta = (p/p0)^(1/equations.gamma)*p0/R
    theta = rho_theta/rho - 300 * exp(0.01^2 / 9.81 * z)
	return theta
end


function compute_global_index(element, ilocal, jlocal, polydeg, cell_x, cell_y)

    i = polydeg * (cell_x-1) + ilocal
    j = polydeg * (cell_y-1) + jlocal

    return i, j
end