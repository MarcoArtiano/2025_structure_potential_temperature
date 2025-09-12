using CairoMakie
using DataFrames
using CSV

function plot_time_momentum(; cells_per_dimension = (100, 60), polydeg = 3, alfa, namenc = "EulerNC", name = "Euler", form = "false", xr_B = 60000)

    t1nc,t2nc,t3nc,t4nc,t5nc,znc = retrieve_data_momentum(namenc, cells_per_dimension, polydeg, alfa, xr_B, form)

    t1,t2,t3,t4,t5,z = retrieve_data_momentum(name, cells_per_dimension, polydeg, alfa, xr_B, form)
    fig = CairoMakie.Figure(size = (1000,600))
    kwargs = (xlabel = L"\overline{m}(z)", xlabelsize = 25,  limits = ((0.5, 1.2), (0, 12)), xticklabelsize = 17.0, yticklabelsize = 17.0,  titlesize = 20)
   
    if name == "Potential"
	ax1 = Axis(fig[2, 1]; kwargs..., title=L"Point-wise $(\varrho, \varrho v, \varrho \theta)$", ylabel = L"$z$ [km]",ylabelsize = 25) 

    ax2 = Axis(fig[2, 2]; kwargs..., title=L"Non-Conservative $(\varrho, \varrho v, \varrho \theta)$")
    elseif name == "Euler"
        ax1 = Axis(fig[2, 1]; kwargs..., title=L"Point-wise $(\varrho, \varrho v, \varrho E)$", ylabel = L"$z$ [km]",ylabelsize = 25) 

    ax2 = Axis(fig[2, 2]; kwargs..., title=L"Non-Conservative $(\varrho, \varrho v, \varrho E)$")
    end
    kwargs_plot = (linewidth = 2.6,)

	lines!(ax1,  t1[:,2], z; label = L"$2.5 h$", kwargs_plot..., color = colors[1])
	lines!(ax1,  t2[:,2], z; label = L"$5 h$", kwargs_plot..., color = colors[2])
	lines!(ax1,  t3[:,2], z; label = L"$7.5 h$", kwargs_plot..., color = colors[3])
	lines!(ax1,  t4[:,2], z; label = L"$10 h$", kwargs_plot..., color = colors[4])
	lines!(ax1,  t5[:,2], z; label = L"$12.5 h$", kwargs_plot..., color = colors[5])

    lines!(ax2,  t1nc[:,2], znc; label = L"$2.5 h$", kwargs_plot..., color = colors[1])
	lines!(ax2,  t2nc[:,2], znc; label = L"$5 h$", kwargs_plot...,color = colors[2])
	lines!(ax2,  t3nc[:,2], znc; label = L"$7.5 h$", kwargs_plot...,color = colors[3])
	lines!(ax2,  t4nc[:,2], znc; label = L"$10 h$", kwargs_plot...,color = colors[4])
	lines!(ax2,  t5nc[:,2], znc; label = L"$12.5 h$", kwargs_plot...,color = colors[5])
    titlelegend = L"\text{Normalized momentum } m = \int_{-\infty}^{\infty}\overline{\varrho}(z)u'(x,z)w'(x,z)dx"
    leg = Legend(fig[1, 1:2], ax1, nothing , titlesize = 24, orientation = :horizontal, labelsize = 25.0)

	save(pwd()*"/test_cases/linearhydrostatic/plots/linearhydrostatic_$(name).pdf", fig)

	return nothing
end


function plot_time_comparison(; cells_per_dimension = (100, 60), polydeg = 3, alfa, form = "false", xr_B = 60000)
    namenc = "EulerNC"
    t1nce,t2nce,t3nce,t4nce,t5nce,znce = retrieve_data_momentum(namenc, cells_per_dimension, polydeg, alfa, xr_B, form)

    namenc = "PotentialNC"
    t1ncp,t2ncp,t3ncp,t4ncp,t5ncp,zncp = retrieve_data_momentum(namenc, cells_per_dimension, polydeg, alfa, xr_B, form)

    name = "Euler"
    t1e,t2e,t3e,t4e,t5e,ze = retrieve_data_momentum(name, cells_per_dimension, polydeg, alfa, xr_B, form)
     name = "Potential"
    t1p,t2p,t3p,t4p,t5p,z = retrieve_data_momentum(name, cells_per_dimension, polydeg, alfa, xr_B, form)
    fig = CairoMakie.Figure(size = (1000,600))
    kwargs = (xlabel = L"\overline{m}(z)", xlabelsize = 25,  limits = ((0.92, 1.08), (0, 12)), xticklabelsize = 17.0, yticklabelsize = 17.0,  titlesize = 20)
   
   
	ax1 = Axis(fig[2, 1]; kwargs..., ylabel = L"$z$ [Km]",ylabelsize = 25) 

    kwargs_plot = (linewidth = 3.5,)

	lines!(ax1,  t5nce[:,2], z; label = L"$\varrho E$ NC", kwargs_plot..., color = colors[1])
    lines!(ax1, t5e[:,2], z; color = colors[2],  label = L"$\varrho E$", kwargs_plot..., linestyle = :solid)
    # scatterlines!(ax1, t5e[:,2], z; marker = :utriangle, color = colors[2], markersize = 10, label = L"$\varrho E$", kwargs_plot...)
    lines!(ax1,  t5ncp[:,2], z; label = L"$\varrho \theta$ NC", kwargs_plot..., color = colors[3], linestyle = :dot)
    # scatterlines!(ax1,  t5p[:,2], z; label = L"$\varrho \theta$",  marker = :circle , markersize = 8, kwargs_plot..., color = colors[4])
    lines!(ax1,  t5p[:,2], z; label = L"$\varrho \theta$", kwargs_plot..., color = colors[4], linestyle = :dash)

    titlelegend = L"\text{Normalized momentum } m = \int_{-\infty}^{\infty}\overline{\varrho}(z)u'(x,z)w'(x,z)dx"
    leg = Legend(fig[1, 1], ax1, nothing , titlesize = 24, orientation = :horizontal, labelsize = 25.0)

	save(pwd()*"/test_cases/linearhydrostatic/plots/linearhydrostatic_comparison.pdf", fig)

	return nothing
end


function plot_time_comparison2(; cells_per_dimension = (100, 60), polydeg = 3, alfa, form = "false", xr_B = 60000)
    namenc = "EulerNC"
    t1nce,t2nce,t3nce,t4nce,t5nce,znce = retrieve_data_momentum(namenc, cells_per_dimension, polydeg, alfa, xr_B, form)

    namenc = "PotentialNC"
    t1ncp,t2ncp,t3ncp,t4ncp,t5ncp,zncp = retrieve_data_momentum(namenc, cells_per_dimension, polydeg, alfa, xr_B, form)

    name = "Euler"
    t1e,t2e,t3e,t4e,t5e,ze = retrieve_data_momentum(name, cells_per_dimension, polydeg, alfa, xr_B, form)
     name = "Potential"
    t1p,t2p,t3p,t4p,t5p,z = retrieve_data_momentum(name, cells_per_dimension, polydeg, alfa, xr_B, form)
    fig = CairoMakie.Figure(size = (1000,600))
    kwargs = (xlabel = L"\overline{m}(z)", xlabelsize = 25,  limits = ((0.9, 1.1), (0, 12)), xticklabelsize = 17.0, yticklabelsize = 17.0,  titlesize = 20)
   
    titlesize = 28

	ax1 = Axis(fig[2, 1]; kwargs..., ylabel = L"$z$ [km]",ylabelsize = 25, title = title=L"$(\varrho, \varrho v, \varrho E)$", titlesize = titlesize) 
    ax2 = Axis(fig[2, 2]; kwargs..., ylabel = L"$z$ [km]",ylabelsize = 25, title = title=L"$(\varrho, \varrho v, \varrho \theta)$", titlesize = titlesize) 

    kwargs_plot = (linewidth = 3.5,)
    lines!(ax1, t5e[:,2], z; color = colors[2],  label = L"Point-wise $$", kwargs_plot..., linestyle = :solid)

	lines!(ax1,  t5nce[:,2], z; label = L"Non-Conservative $$", kwargs_plot..., color = colors[1])
    # scatterlines!(ax1, t5e[:,2], z; marker = :utriangle, color = colors[2], markersize = 10, label = L"$\varrho E$", kwargs_plot...)
       lines!(ax2,  t5p[:,2], z; label = L"Point-wise $$", kwargs_plot..., color = colors[2], linestyle = :solid)

    lines!(ax2,  t5ncp[:,2], z; label = L"Non-Conservative $$", kwargs_plot..., color = colors[1], linestyle = :solid)
    # scatterlines!(ax1,  t5p[:,2], z; label = L"$\varrho \theta$",  marker = :circle , markersize = 8, kwargs_plot..., color = colors[4])
    titlelegend = L"\text{Normalized momentum } m = \int_{-\infty}^{\infty}\overline{\varrho}(z)u'(x,z)w'(x,z)dx"
    leg = Legend(fig[1, 1:2], ax1, nothing , titlesize = 27, orientation = :horizontal, labelsize = 25.0)

	save(pwd()*"/test_cases/linearhydrostatic/plots/linearhydrostatic_comparison_2.pdf", fig)

	return nothing
end


function retrieve_data_momentum(name, cells_per_dimension, polydeg, alfa, xr_B, form)
    folder_path = pwd()*"/test_cases/linearhydrostatic/out/"

    t1 = CSV.read(folder_path*"$(name)_2.5_$(cells_per_dimension[1])x$(cells_per_dimension[2])_$(polydeg)_$(alfa)_$(xr_B)_$(form).csv", DataFrame)
	t2 = CSV.read(folder_path*"$(name)_5_$(cells_per_dimension[1])x$(cells_per_dimension[2])_$(polydeg)_$(alfa)_$(xr_B)_$(form).csv", DataFrame)
	t3 = CSV.read(folder_path*"$(name)_7.5_$(cells_per_dimension[1])x$(cells_per_dimension[2])_$(polydeg)_$(alfa)_$(xr_B)_$(form).csv", DataFrame)
	t4 = CSV.read(folder_path*"$(name)_10_$(cells_per_dimension[1])x$(cells_per_dimension[2])_$(polydeg)_$(alfa)_$(xr_B)_$(form).csv", DataFrame)
	t5 = CSV.read(folder_path*"$(name)_12.5_$(cells_per_dimension[1])x$(cells_per_dimension[2])_$(polydeg)_$(alfa)_$(xr_B)_$(form).csv", DataFrame)
	z = t1[:,1]./1e3
    return t1,t2,t3,t4,t5,z
end

colors = Makie.wong_colors()

plot_time_momentum(alfa = 0.035)

plot_time_momentum(namenc = "PotentialNC", name = "Potential", alfa = 0.035)

plot_time_comparison(alfa=0.035)
plot_time_comparison2(alfa=0.035)