using CSV, DataFrames, CairoMakie, LaTeXStrings, PotentialTemperature

function convergence_plots_warped!(vars, fluxes, yboundss, title, filename)
     flux_noncons = flux_nonconservative_gravity_log
     var = vars[1]
     flux = fluxes[1]
     ybounds = yboundss[1]
    fig = CairoMakie.Figure(size=(1600, 800))
    CFL = 0.1
    df1 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)Euler_flux_ranocha___.csv", DataFrame)
    df2 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)EUNC_zero_flux_ranocha_$(flux_noncons)_.csv", DataFrame)

    df3 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)PT_" * flux * "___.csv", DataFrame)
    df4 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)PT_NC_zero_" * flux * "_$(flux_noncons)_.csv", DataFrame)

    ymin = ybounds[1]
    ymax = ybounds[2]
    xlabel1 = L"$\Delta x$ [km]"
        if var == "wp"
    ylabel1 = L"||w||_2"
    titlelegend = L"L_2 \text{ Error } w"
    elseif var == "up"
        ylabel1 = L"||u||_2"
        titlelegend = L"L_2 \text{ Error } u"
    elseif var == "pp"
        ylabel1 = L"||p||_2"
        titlelegend = L"L_2 \text{ Error } p"
    elseif var == "Tp"
        ylabel1 = L"||T||_2"
        titlelegend = L"L_2 \text{ Error } T"
    elseif var == "thetap"
        ylabel1 = L"||\theta||_2"
        titlelegend = L"L_2 \text{ Error } \theta"
    end
    titlesize = 27
    xlabelsize = 22.0
    ylabelsize = 27
    xticklabelsize = 19.0
    yticklabelsize = 19.0
    xticks1 = ([1.25, 2.5, 5.0, 10.0])
    yticks1 = LogTicks(WilkinsonTicks(6, k_min=4))
  
    ax1 = Axis(fig[2, 1], title=L"Point-wise $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax2 = Axis(fig[2, 2], title=L"Non-Conservative $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax3 = Axis(fig[3, 1], title=L"Point-wise $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax4 = Axis(fig[3, 2], title=L"Non-Conservative $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)

    plot_data!(ax1, df1)
    plot_data!(ax2, df2)
    plot_data!(ax3, df3)
    plot_data!(ax4, df4)
    leg = Legend(fig[1, 1:2], ax3, titlelegend, titlesize = 30, orientation = :horizontal, labelsize = 28.0)
    var = vars[2]
     flux = fluxes[2]
     ybounds = yboundss[2]
     df1 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)Euler_flux_ranocha___.csv", DataFrame)
    df2 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)EUNC_zero_flux_ranocha_$(flux_noncons)_.csv", DataFrame)

    df3 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)PT_" * flux * "___.csv", DataFrame)
    df4 = CSV.read(pwd() * "/test_cases/gravitywaves/out/warped_" * var * "_error_CFL_$(CFL)PT_NC_zero_" * flux * "_$(flux_noncons)_.csv", DataFrame)

    ymin = ybounds[1]
    ymax = ybounds[2]
    xlabel1 = L"$\Delta x$ [km]"
        if var == "wp"
    ylabel1 = L"||w||_2"
    titlelegend = L"L_2 \text{ Error } w"
    elseif var == "up"
        ylabel1 = L"||u||_2"
        titlelegend = L"L_2 \text{ Error } u"
    elseif var == "pp"
        ylabel1 = L"||p||_2"
        titlelegend = L"L_2 \text{ Error } p"
    elseif var == "Tp"
        ylabel1 = L"||T||_2"
        titlelegend = L"L_2 \text{ Error } T"
    elseif var == "thetap"
        ylabel1 = L"||\theta||_2"
        titlelegend = L"L_2 \text{ Error } \theta"
    end
    xlabelsize = 20.0
    ylabelsize = 25
    xticklabelsize = 17.0
    yticklabelsize = 17.0
    xticks1 = ([1.25, 2.5, 5.0, 10.0])
    yticks1 = LogTicks(WilkinsonTicks(6, k_min=4))
  
    ax1 = Axis(fig[2, 3], title=L"Point-wise $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax2 = Axis(fig[2, 4], title=L"Non-Conservative $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax3 = Axis(fig[3, 3], title=L"Point-wise $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax4 = Axis(fig[3, 4], title=L"Non-Conservative $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)

    plot_data!(ax1, df1)
    plot_data!(ax2, df2)
    plot_data!(ax3, df3)
    plot_data!(ax4, df4)
    leg = Legend(fig[1, 3:4], ax3, titlelegend, titlesize = 30, orientation = :horizontal, labelsize = 28.0)

    save(pwd()*"/test_cases/gravitywaves/plots/warped_"*filename*"_"*vars[1]* vars[2]*".pdf", fig)
end

function convergence_plots!(vars, fluxes, ybounds, title, filename; CFL = 0.1, flux_noncons = flux_nonconservative_gravity_log)
    var1 = vars[1]
    var2 = vars[2]

    flux1 = fluxes[1]
    flux2 = fluxes[2]

    ybounds1 = ybounds[1]
    ybounds2 = ybounds[2]
    titlesizeleg = 30
    titlesize = 27
    fig = CairoMakie.Figure(size=(1600, 800))

    df1 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var1 * "_error_CFL_$(CFL)Euler_flux_ranocha___.csv", DataFrame)
    df2 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var1 * "_error_CFL_$(CFL)EUNC_flux_ranocha_$(flux_noncons)_.csv", DataFrame)

    df3 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var1 * "_error_CFL_$(CFL)PT_" * flux1 * "___.csv", DataFrame)
    df4 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var1 * "_error_CFL_$(CFL)PT_NC_" * flux1 * "_$(flux_noncons)_.csv", DataFrame)

    ymin = ybounds1[1]
    ymax = ybounds1[2]
    xlabel1 = L"$\Delta x$ [km]"
    if var1 == "wp"
    ylabel1 = L"||w||_2"
    titlelegend = L"L_2 \text{ Error } w"
    elseif var1 == "up"
        ylabel1 = L"||u||_2"
        titlelegend = L"L_2 \text{ Error } u"
    elseif var1 == "pp"
        ylabel1 = L"||p||_2"
        titlelegend = L"L_2 \text{ Error } p"
    elseif var1 == "Tp"
        ylabel1 = L"||T||_2"
        titlelegend = L"L_2 \text{ Error } T"
    elseif var1 == "thetap"
        ylabel1 = L"||\theta||_2"
        titlelegend = L"L_2 \text{ Error } \theta"
    end
    xlabelsize = 22.0
    ylabelsize = 27
    xticklabelsize = 19.0
    yticklabelsize = 19.0

    xticks1 = ([1.25, 2.5, 5.0, 10.0])
    yticks1 = LogTicks(WilkinsonTicks(6, k_min=4))
    ax1 = Axis(fig[2, 1], title=L"Point-wise $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax2 = Axis(fig[2, 2], title=L"Non-Conservative $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax3 = Axis(fig[3, 1], title=L"Point-wise $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax4 = Axis(fig[3, 2], title=L"Non-Conservative $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)

    plot_data!(ax1, df1)
    plot_data!(ax2, df2)
    plot_data!(ax3, df3)
    plot_data!(ax4, df4)
    leg = Legend(fig[1, 1:2], ax3, titlelegend, titlesize = titlesizeleg, orientation = :horizontal, labelsize = 28.0)

    df1 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var2 * "_error_CFL_$(CFL)Euler_flux_ranocha___.csv", DataFrame)
    df2 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var2 * "_error_CFL_$(CFL)EUNC_flux_ranocha_$(flux_noncons)_.csv", DataFrame)

    df3 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var2 * "_error_CFL_$(CFL)PT_" * flux2 * "___.csv", DataFrame)
    df4 = CSV.read(pwd() * "/test_cases/gravitywaves/out/" * var2 * "_error_CFL_$(CFL)PT_NC_" * flux2 * "_$(flux_noncons)_.csv", DataFrame)

    ymin = ybounds2[1]
    ymax = ybounds2[2]
    xlabel1 = L"$\Delta x$ [km]"
    if var2 == "wp"
    ylabel1 = L"||w||_2"
    titlelegend = L"L_2 \text{ Error } w"
    elseif var2 == "up"
        ylabel1 = L"||u||_2"
        titlelegend = L"L_2 \text{ Error } u"
    elseif var2 == "pp"
        ylabel1 = L"||p||_2"
        titlelegend = L"L_2 \text{ Error } p"
    elseif var2 == "Tp"
        ylabel1 = L"||T||_2"
        titlelegend = L"L_2 \text{ Error } T"
    elseif var2 == "thetap"
        ylabel1 = L"||\theta||_2"
        titlelegend = L"L_2 \text{ Error } \theta"
    end
    xticks1 = ([1.25, 2.5, 5.0, 10.0])
    yticks1 = LogTicks(WilkinsonTicks(6, k_min=4))
    kwargs = (limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, ylabel = ylabel1, xscale=log10, yscale=log10, xlabelsize = 30, ylabelsize = 25, xticklabelsize = 17.0, yticklabelsize = 17.0)
    ax1 = Axis(fig[2, 3], title=L"Point-wise $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax2 = Axis(fig[2, 4], title=L"Non-Conservative $(\varrho, \varrho v, \varrho E)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax3 = Axis(fig[3, 3], title=L"Point-wise $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1,  xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)
    ax4 = Axis(fig[3, 4], title=L"Non-Conservative $(\varrho, \varrho v, \varrho \theta)$", limits=(nothing, (ymin, ymax)), xticks=xticks1, yticks=yticks1, xlabel = xlabel1, xscale=log10, yscale=log10, xlabelsize = xlabelsize, ylabelsize = ylabelsize, xticklabelsize = xticklabelsize, yticklabelsize = yticklabelsize, titlesize = titlesize)

    plot_data!(ax1, df1)
    plot_data!(ax2, df2)
    plot_data!(ax3, df3)
    plot_data!(ax4, df4)
    leg = Legend(fig[1, 3:4], ax3, titlelegend, titlesize = titlesizeleg, orientation = :horizontal, labelsize = 28.0)

    save(pwd()*"/test_cases/gravitywaves/plots/"*filename*"_"*var1*var2*".pdf", fig)
end

function plot_data!(ax, df)

    unique_degrees = sort(unique(df.polydeg))
    for (i, deg) in enumerate(unique_degrees)
        subdf = filter(row -> row.polydeg == deg, df)
        scatterlines!(ax, subdf.dx ./ 1000, subdf.l2_error,
            color=colors[i], label=L"$p = %$(deg)$", markersize=17,  marker = markers[i])
        x_ref = [maximum(df.dx), minimum(df.dx)]
        if deg == 4 #&& df.problem_name[1] == "EUNC"
            y_ref = subdf.l2_error[1] * ((x_ref) ./ subdf.dx[1]) .^ (deg + 1)
        else
            y_ref = subdf.l2_error[end] * ((x_ref) ./ subdf.dx[end]) .^ (deg + 1)
        end
        order = deg + 1
        label = latexstring("\$\\mathcal{O}(\\Delta x^{ $order })\$")  # con i $ per inline math
                #label = L"Order  %$(order)$"
        lines!(ax, x_ref ./ 1000, y_ref, color=colors[i], linestyle=linestyles[i], label=label, linewidth = 3)
    end
end

colors = Makie.wong_colors()
markers = (:circle, :star5, :rect)
linestyles = (:dash, :dashdot, :dashdotdot)
convergence_plots!(("up","wp"), ("flux_theta","flux_theta"), ((10^-10, 10^-3), (10^-10, 10^0)), "L2 error w", "Convergencefluxtheta")
convergence_plots!(("pp","Tp"), ("flux_theta","flux_theta"), ((10^-10, 10^3), (10^-10, 10^1)), "L2 error w", "Convergencefluxtheta")

convergence_plots_warped!(("wp","Tp"), ("flux_theta","flux_theta"), ((10^-10, 10^1),(10^-10, 10^2)), "Convergence of w", "Convergencefluxtheta")