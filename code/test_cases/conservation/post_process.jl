using DelimitedFiles
using CairoMakie
using LaTeXStrings

function load_processed_data(file_name, folder_path)
	file_path = joinpath(folder_path, file_name)
	return readdlm(file_path)
end

function tgv_plots(initial_refinement_level, polydeg)

    folder_path = pwd()*"/test_cases/conservation/out/tgv"
	euler = load_processed_data("euler_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta = load_processed_data("potential_flux_theta_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_am = load_processed_data("potential_flux_theta_AM_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_rhos = load_processed_data("potential_flux_theta_rhos_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_rhos_am = load_processed_data("potential_flux_theta_rhos_AM_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_global = load_processed_data("potential_flux_theta_global_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)

    data_tuple = (euler, theta, theta_am, theta_rhos, theta_rhos_am, theta_global)
	label_tuple = (L"\text{Ranocha}", L"\text{TEC log}", L"\text{TEC avg}", L"\text{EC log}", L"\text{EC avg}", L"\text{ETEC}")
    linestyle_tuple = (:solid, :dashdotdot, :dash, :dash, :dot, :dashdot)

    fig = Figure(size = (1000,500))
    kwargs = (xlabel = L"t", xlabelsize = 30, ylabelsize = 25, limits = ((0, 50), nothing), xticklabelsize = 17.0, yticklabelsize = 17.0)
    ylabel = L"\langle \varrho s \rangle"
	ax1 = Axis(fig[2, 1]; ylabel = ylabel, kwargs...)

    ylabel = L"\langle \varrho E \rangle"
    ax2 = Axis(fig[2, 2]; ylabel = ylabel, kwargs... )
    index = 4
    create_plot!(data_tuple, label_tuple,linestyle_tuple, ylabel, index, ax1, position = :lb)

    index = 3
    create_plot!(data_tuple, label_tuple, linestyle_tuple, ylabel, index, ax2, position = :lt)
    Label(fig[2, 1, Top()], halign = :left, "×10⁻⁷", fontsize = 15.0)
    Label(fig[2, 2, Top()], halign = :left, "×10⁻⁷", fontsize = 15.0)
    leg = Legend(fig[1, 1:2], ax1, orientation = :horizontal, labelsize = 20.0)
    save(pwd()*"/test_cases/conservation/plots/tgv.pdf", fig)
    save(pwd()*"/test_cases/conservation/plots/tgv.png", fig)

    return nothing
end

function dw_plots(initial_refinement_level, polydeg)

    folder_path = pwd()*"/test_cases/conservation/out/density_wave"
	euler = load_processed_data("euler_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta = load_processed_data("potential_flux_theta_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_am = load_processed_data("potential_flux_theta_AM_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_rhos = load_processed_data("potential_flux_theta_rhos_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_rhos_am = load_processed_data("potential_flux_theta_rhos_AM_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)
	theta_global = load_processed_data("potential_flux_theta_global_ref$(initial_refinement_level)_polydeg$polydeg.dat", folder_path)

    data_tuple = (euler, theta, theta_am, theta_rhos, theta_rhos_am, theta_global)
	label_tuple = (L"\text{Ranocha}", L"\text{TEC log}", L"\text{TEC avg}", L"\text{EC log}", L"\text{EC avg}", L"\text{ETEC}")
    linestyle_tuple = (:dash, :solid, :solid, :solid, :solid, :solid)
	linestyle_tuple = (:solid, :dashdotdot, :dash, :dash, :dot, :dashdot)

    fig = Figure(size = (1000,500))
    kwargs = (xlabel = L"t", xlabelsize = 30, ylabelsize = 25, limits = ((0, 40), nothing), xticklabelsize = 17.0, yticklabelsize = 17.0)
    ylabel = L"\langle \varrho s \rangle"
	ax1 = Axis(fig[2, 1]; ylabel = ylabel, kwargs...)

    ylabel = L"\langle \varrho E \rangle"
    ax2 = Axis(fig[2, 2]; ylabel = ylabel, kwargs... )
    index = 4
    create_plot!(data_tuple, label_tuple,linestyle_tuple, ylabel, index, ax1, position = :lb)

    index = 3
    create_plot!(data_tuple, label_tuple, linestyle_tuple, ylabel, index, ax2, position = :lt)
    Label(fig[2, 1, Top()], halign = :left, "×10⁻⁷", fontsize = 15.0)
    Label(fig[2, 2, Top()], halign = :left, "×10⁻⁷", fontsize = 15.0)
    leg = Legend(fig[1, 1:2], ax1, orientation = :horizontal, labelsize = 20.0)
    save(pwd()*"/test_cases/conservation/plots/dw.pdf", fig)
    save(pwd()*"/test_cases/conservation/plots/dw.png", fig)

    return nothing
end

function create_plot!(data_tuple, label_tuple, linestyle_tuple, ylabel, index, ax; position = :lt)

	for (var, label, linestyle) in zip(data_tuple, label_tuple, linestyle_tuple)
		plot_var!(var, index, label, ax, linestyle)
	end
	return nothing
end

function plot_var!(var, index, legend, ax, linestyle)

	lines!(ax, var[:,1], var[:, index]./10^(-7); label = legend, linewidth = 3.5, linestyle = linestyle)

end

tgv_plots(5,0)
dw_plots(6,0)