using PyPlot

matplotlib = PyPlot.matplotlib
rcParams = matplotlib["rcParams"]

rcParams["pdf.fonttype"] = 42  # Usa Type 42 fonts TrueType, che vengono embedded subset nel PDF
rcParams["ps.fonttype"] = 42

rcParams["text.usetex"] = true  # Se vuoi usare LaTeX per i testi, altrimenti false
rcParams["font.family"] = "serif"  #

@inline function contour_baroclinic(sol_euler, semi_euler,sol_theta, semi_theta, (Kh, Kv), nvisnodes, equations_euler, equations_theta)
    
    data_layer_euler, nodes_layer_euler = retrieve_values_per_layer(Kh, Kv, semi_euler, sol_euler)
    data_layer_theta, nodes_layer_theta = retrieve_values_per_layer(Kh, Kv, semi_theta, sol_theta)

	p_euler = retrieve_pressure_layer(data_layer_euler, equations_euler)
	T_euler = retrieve_temperature_layer(data_layer_euler, equations_euler)

    p_theta = retrieve_pressure_layer(data_layer_theta, equations_theta)
	T_theta = retrieve_temperature_layer(data_layer_theta, equations_theta)

	p_euler, x, y, z = spectral_interpolation3D(p_euler, nodes_layer_euler, semi_euler, nvisnodes)
    p_theta, x, y, z = spectral_interpolation3D(p_theta, nodes_layer_theta, semi_theta, nvisnodes)

	T_euler, _, _, _ = spectral_interpolation3D(T_euler, nodes_layer_euler, semi_euler, nvisnodes)
    T_theta, _, _, _ = spectral_interpolation3D(T_theta, nodes_layer_theta, semi_theta, nvisnodes)

	lon, lat = cart2sphere(x, y, z)

	lonshift = 60
	lon = @. mod(lon + 180 - lonshift, 360) - 180

	mask = (lat .>= 0) .& (lon .>= -lonshift)

	lon = lon[mask]
	lat = lat[mask]
	p_euler = p_euler[mask]
	p_theta = p_theta[mask]
	t_euler = T_euler[mask]
	t_theta = T_theta[mask]
	plevels = vcat([955], [960 + 5i for i ∈ 0:11], [1020])
	plevels = 11
	pnorm = matplotlib.colors.TwoSlopeNorm(vmin = 955, vcenter = 990, vmax = 1025)
	pnorm = nothing
	cmap = ColorMap("plasma").copy()
	shrinkcb = 0.7
    fig, axs = subplots(1, 1, figsize=(27, 20))
	cs1 = axs.tricontour(lon, lat, p_euler; levels = plevels, colors = ("k",))
	cset = axs.tricontourf(
		lon,
		lat,
		p_theta;
		levels = plevels,
		cmap,
		norm = pnorm,
		extend = "neither",
	)

    axs[:set_title](raw"Day 10 $(\varrho, \varrho v, \varrho E)$, Surface Pressure [hPa]", fontsize=40)

    cbar_ax = fig.add_axes([0.94, 0.25, 0.015, 0.5])  # Cambia qui per regolare posizione
    cbar = colorbar(
		cset,
        cax = cbar_ax,
        ax = axs,
		ticks = plevels isa Int ? nothing : plevels[1+2:2:end-1],
		shrink = shrinkcb)
    plt.subplots_adjust(wspace = 0.05)
    axs.tick_params(axis="both",labelsize=40)  # Sostituisci 14 con la dimensione desiderata

  xticks = [-60, -30, 0, 30, 60, 90, 120, 150, 180]
  xticklabels =
    ["0", "30E", "60E", "90E", "120E", "150E", "180", "150W", "120W"]
  yticks = [0, 30, 60, 90]
  yticklabels = ["0", "30N", "60N", "90N"]
    axs.set_xticks(xticks)
    axs.set_xticklabels(xticklabels)
    axs.set_yticks(yticks)
    axs.set_yticklabels(yticklabels)
    axs.set_xlim([xticks[1], xticks[end]])
    axs.set_ylim([yticks[1], yticks[end]])
    axs.set_aspect(1)
  cbar[:ax][:tick_params](labelsize=30)

	savefig(joinpath("test_cases/baroclinic/plots/contour_pressure_euler_new_colormap_$(Kh)x$(Kv).png"),dpi = 600, bbox_inches="tight")
	close(fig)
	
	    fig, axs = subplots(1, 1, figsize=(27, 20))
	cs1 = axs.tricontour(lon, lat, t_euler; levels = plevels, colors = ("k",))
	cset = axs.tricontourf(
		lon,
		lat,
		t_euler;
		levels = plevels,
		cmap,
		norm = pnorm,
		extend = "neither",
	)

    axs[:set_title](raw"Day 10 $(\varrho, \varrho v, \varrho E)$, Temperature [K]", fontsize=40)

    cbar_ax = fig.add_axes([0.94, 0.25, 0.015, 0.5])  # Cambia qui per regolare posizione
    cbar = colorbar(
		cset,
        cax = cbar_ax,
        ax = axs,
		ticks = plevels isa Int ? nothing : plevels[1+2:2:end-1],
		shrink = shrinkcb)
    plt.subplots_adjust(wspace = 0.05)
    axs.tick_params(axis="both",labelsize=40)  # Sostituisci 14 con la dimensione desiderata

  xticks = [-60, -30, 0, 30, 60, 90, 120, 150, 180]
  xticklabels =
    ["0", "30E", "60E", "90E", "120E", "150E", "180", "150W", "120W"]
  yticks = [0, 30, 60, 90]
  yticklabels = ["0", "30N", "60N", "90N"]
    axs.set_xticks(xticks)
    axs.set_xticklabels(xticklabels)
    axs.set_yticks(yticks)
    axs.set_yticklabels(yticklabels)
    axs.set_xlim([xticks[1], xticks[end]])
    axs.set_ylim([yticks[1], yticks[end]])
    axs.set_aspect(1)
  cbar[:ax][:tick_params](labelsize=30)

	savefig(joinpath("test_cases/baroclinic/plots/contour_temperature_euler_new_colormap_$(Kh)x$(Kv).png"),dpi = 600, bbox_inches="tight")
	close(fig)


    fig, axs = subplots(1, 1, figsize=(27, 20))
	cs1 = axs.tricontour(lon, lat, p_theta; levels = plevels, colors = ("k",))
	cset = axs.tricontourf(
		lon,
		lat,
		p_theta;
		levels = plevels,
		cmap,
		norm = pnorm,
		extend = "neither",
	)

    axs[:set_title](raw"Day 10 $(\varrho, \varrho v, \varrho \theta)$, Surface Pressure [hPa]", fontsize=40)

    cbar_ax = fig.add_axes([0.94, 0.25, 0.015, 0.5])
    cbar = colorbar(
		cset,
        cax = cbar_ax,
        ax = axs,
		ticks = plevels isa Int ? nothing : plevels[1+2:2:end-1],
		shrink = shrinkcb)
    plt.subplots_adjust(wspace = 0.05)
    axs.tick_params(axis="both",labelsize=40)

  xticks = [-60, -30, 0, 30, 60, 90, 120, 150, 180]
  xticklabels =
    ["0", "30E", "60E", "90E", "120E", "150E", "180", "150W", "120W"]
  yticks = [0, 30, 60, 90]
  yticklabels = ["0", "30N", "60N", "90N"]
    axs.set_xticks(xticks)
    axs.set_xticklabels(xticklabels)
    axs.set_yticks(yticks)
    axs.set_yticklabels(yticklabels)
    axs.set_xlim([xticks[1], xticks[end]])
    axs.set_ylim([yticks[1], yticks[end]])
    axs.set_aspect(1)
  cbar[:ax][:tick_params](labelsize=30)

	savefig(joinpath("test_cases/baroclinic/plots/contour_pressure_theta_new_colormap_$(Kh)x$(Kv).png"),dpi = 600, bbox_inches="tight")
	close(fig)

		    fig, axs = subplots(1, 1, figsize=(27, 20))
	cs1 = axs.tricontour(lon, lat, t_theta; levels = plevels, colors = ("k",))
	cset = axs.tricontourf(
		lon,
		lat,
		t_theta;
		levels = plevels,
		cmap,
		norm = pnorm,
		extend = "neither",
	)

    axs[:set_title](raw"Day 10 $(\varrho, \varrho v, \varrho \theta)$, Temperature [K]", fontsize=40)

    cbar_ax = fig.add_axes([0.94, 0.25, 0.015, 0.5])  # Cambia qui per regolare posizione
    cbar = colorbar(
		cset,
        cax = cbar_ax,
        ax = axs,
		ticks = plevels isa Int ? nothing : plevels[1+2:2:end-1],
		shrink = shrinkcb)
    plt.subplots_adjust(wspace = 0.05)
    axs.tick_params(axis="both",labelsize=40)  # Sostituisci 14 con la dimensione desiderata

  xticks = [-60, -30, 0, 30, 60, 90, 120, 150, 180]
  xticklabels =
    ["0", "30E", "60E", "90E", "120E", "150E", "180", "150W", "120W"]
  yticks = [0, 30, 60, 90]
  yticklabels = ["0", "30N", "60N", "90N"]
    axs.set_xticks(xticks)
    axs.set_xticklabels(xticklabels)
    axs.set_yticks(yticks)
    axs.set_yticklabels(yticklabels)
    axs.set_xlim([xticks[1], xticks[end]])
    axs.set_ylim([yticks[1], yticks[end]])
    axs.set_aspect(1)
  cbar[:ax][:tick_params](labelsize=30)

	savefig(joinpath("test_cases/baroclinic/plots/contour_temperature_theta_new_colormap_$(Kh)x$(Kv).png"),dpi = 600, bbox_inches="tight")
	close(fig)

end

function retrieve_pressure_layer(data_layer, equations)
	size_ = size(data_layer)
	p = zeros(Float64, size_[2:end]...)
	for i in 1:size_[2]
		for j in 1:size_[3]
			for element in 1:size_[4]
				p[i, j, element] = cons2pressure(data_layer[:, i, j, element], equations)
			end
		end
	end
	return p

end

function retrieve_temperature_layer(data_layer, equations)
	size_ = size(data_layer)
	T = zeros(Float64, size_[2:end]...)
	for i in 1:size_[2]
		for j in 1:size_[3]
			for element in 1:size_[4]
				T[i, j, element] = cons2temperature(data_layer[:, i, j, element], equations)
			end
		end
	end
	return T

end

function cons2temperature(u, equations::Union{CompressibleEulerEquations3D, CompressibleEulerEquations3DNC})
	rho, rho_v1, rho_v2, rho_v3, rho_e = u

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v3 = rho_v3 / rho
	p = (1004 / 717 - 1) * (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

	return p / (rho * 287)
end

function cons2temperature(u, equations::Union{CompressibleEulerPotentialTemperatureEquations3D, CompressibleEulerPotentialTemperatureEquations3DNC})
	rho, rho_v1, rho_v2, rho_v3, rho_theta = u

	p = equations.K * rho_theta^equations.gamma
	return p / (rho * 287)
end

function cons2pressure(u, equations::Union{CompressibleEulerEquations3D, CompressibleEulerEquations3DNC})
	rho, rho_v1, rho_v2, rho_v3, rho_e = u

	v1 = rho_v1 / rho
	v2 = rho_v2 / rho
	v3 = rho_v3 / rho
	p = (1004 / 717 - 1) * (rho_e - 0.5f0 * (rho_v1 * v1 + rho_v2 * v2 + rho_v3 * v3))

	return p / 1e2
end

function cons2pressure(u, equations::Union{CompressibleEulerPotentialTemperatureEquations3D, CompressibleEulerPotentialTemperatureEquations3DNC})
	rho, rho_v1, rho_v2, rho_v3, rho_theta = u

	p = equations.K * rho_theta^equations.gamma
	return p / 1e2
end


function cart2sphere(x, y, z)
	r = similar(x)
	r .= sqrt.(x .^ 2 + y .^ 2 + z .^ 2)

	lat = similar(x)
	lon = similar(x)

	lat = asin.(z ./ r) .* (180 / π)
	lon = atan.(y, x) .* (180 / π)

	return vec(lon), vec(lat)
end

@inline function retrieve_values_per_layer(Kh, Kv, semi, sol; layer = 1)

	node_coordinates = semi.cache.elements.node_coordinates
	data = Trixi.wrap_array(sol.u[end], semi)
	size_sol = size(data)
	size_coords = size(node_coordinates)
	sphere_layer = (Kh^2*(layer-1)+1):Kh^2*Kv:6*Kv*Kh^2
	data_layer = zeros(Float64, size_sol[1:end-2]..., Kh^2 * 6)
	nodes_layer = zeros(Float64, size_coords[1:end-2]..., Kh^2 * 6)

	for block in 1:6
		data_layer[:, :, :, (Kh^2*(block-1)+1):Kh^2*block] .= data[:, :, :, 1, sphere_layer[block]:(sphere_layer[block]+Kh^2-1)]
		nodes_layer[:, :, :, (Kh^2*(block-1)+1):Kh^2*block] .= node_coordinates[:, :, :, 1, sphere_layer[block]:(sphere_layer[block]+Kh^2-1)]
	end
	return data_layer, nodes_layer
end

function spectral_interpolation3D(data_layer, nodes_layer, semi, nvisnodes)
	nvars = 1
	size_ = size(nodes_layer)
	n_nodes_2d = Trixi.nnodes(semi.solver)^2
	n_elements = size_[end]
	#  plotting_interp_matrix = plotting_interpolation_matrix_no_boundary(semi.solver; nvisnodes=nvisnodes)
	plotting_interp_matrix = Trixi.plotting_interpolation_matrix(semi.solver; nvisnodes = nvisnodes)

	uEltype = eltype(data_layer)
	x = reshape(view(nodes_layer, 1, :, :, :), n_nodes_2d,
		n_elements)
	y = reshape(view(nodes_layer, 2, :, :, :), n_nodes_2d,
		n_elements)
	z = reshape(view(nodes_layer, 3, :, :, :), n_nodes_2d,
		n_elements)

	u_extracted = StructArray{SVector{nvars, uEltype}}(ntuple(_ -> similar(x,
			(n_nodes_2d,
				n_elements)),
		nvars))
	for element in 1:n_elements
		sk = 1
		for j in eachnode(semi.solver), i in eachnode(semi.solver)
			u_node = SVector(data_layer[i, j, element])
			u_extracted[sk, element] = u_node
			sk += 1
		end
	end
	uplot = StructArray{SVector{nvars, uEltype}}(map(x -> plotting_interp_matrix * x,
		StructArrays.components(u_extracted)))
	xplot, yplot, zplot = plotting_interp_matrix * x, plotting_interp_matrix * y, plotting_interp_matrix * z

	uplot = StructArray{SVector{nvars, uEltype}}(map(x -> plotting_interp_matrix * x,
		StructArrays.components(u_extracted)))

	return getindex.(vec(uplot), 1), vec(xplot), vec(yplot), vec(zplot)
end

function plotting_interpolation_matrix_no_boundary(dg::DGSEM;
	nvisnodes = 2 * length(dg.basis.nodes))

	dξ = 2 / nvisnodes
	interia = [-1 + (j - 1 / 2) * dξ for j ∈ 1:nvisnodes]
	Vp1D = Trixi.polynomial_interpolation_matrix(dg.basis.nodes, interia)
	# For quadrilateral elements, interpolation to plotting nodes involves applying a 1D interpolation
	# operator to each line of nodes. This is equivalent to multiplying the vector containing all node
	# node coordinates on an element by a Kronecker product of the 1D interpolation operator (e.g., a
	# multi-dimensional interpolation operator).
	return Trixi.kron(Vp1D, Vp1D)
end

