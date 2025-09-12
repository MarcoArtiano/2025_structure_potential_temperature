using DelimitedFiles
using CairoMakie
function load_processed_data(file_name, folder_path)
	file_path = joinpath(folder_path, file_name)
 	return readdlm(file_path)
end
folder_path = pwd()*"/test_cases/balance/out"
dd = load_processed_data("well_balancing_Double64.dat", folder_path)
ff = load_processed_data("well_balancing_Float64.dat", folder_path)
ff_theta = load_processed_data("well_balancing_theta_Float64.dat", folder_path)
dd_theta = load_processed_data("well_balancing_theta_Double64.dat", folder_path)
colors = Makie.wong_colors()
stepsize = 2500
stepmarker = 40000
fig = Figure(size = (1200, 600))
kwargs = (xlabel = L"t", xlabelsize = 24, xticklabelsize = 17.0, yticklabelsize = 17.0,  titlesize = 20)
ax1 = Axis(fig[2, 1]; kwargs..., ylabelsize = 25, title = L"\text{Double64}", limits = ((0 ,5000) ,(-0.25, 4) ))
n = 500000+2
markers = (:circle,)
lines!(ax1, dd[2:n,2], dd[2:n,end-1]./1e-18; label = L"$||u||_2 (T = \text{const})$", linewidth = 2.0, color = colors[1])
lines!(ax1, dd[2:n,2], dd[2:n,end]./1e-18; label = L"$||w||_2 (T = \text{const})$", linewidth = 4.0, color = colors[1], linestyle =:dash)
lines!(ax1, dd_theta[2:n,1], dd_theta[2:n,end-1]./1e-18; label = L"$||u||_2 (\theta = \text{const})$", linewidth = 2.0, color = colors[2])
lines!(ax1, dd_theta[2:n,1], dd_theta[2:n,end]./1e-18; label = L"$||w||_2 (\theta = \text{const})$", linewidth = 4.0, linestyle =:dash, color = colors[2],)

ax2 = Axis(fig[2, 2]; kwargs..., ylabelsize = 25, title = L"\text{Float64}", limits = ((0 ,5000) ,(-0.25, 6) ))
lines!(ax2, ff[2:n,2], ff[2:n,end-1]./1e-8; label = L"$||u||_2$", linewidth = 2.0, color = colors[1])
lines!(ax2, ff[2:n,2], ff[2:n,end]./1e-8; label = L"$||w||_2$", linewidth = 4.0, color = colors[1], linestyle =:dash)
lines!(ax2, ff_theta[2:n,2], ff_theta[2:n,end-1]./1e-8; label = L"$||u||_2$", linewidth = 2.0, color = colors[2])
lines!(ax2, ff_theta[2:n,2], ff_theta[2:n,end]./1e-8; label = L"$||w||_2$", linewidth = 4.0, color = colors[2], linestyle =:dash)

handles_colors = [
    LineElement(color = colors[1], linewidth = 3),
    LineElement(color = colors[2], linewidth = 3), LineElement(color = :black, linewidth = 3, linestyle = :solid),
    LineElement(color = :black, linewidth = 3, linestyle = :dash)
]
labels_colors = [L"T = \text{const}", L"\theta = \text{const}", L"||u||_2", L"||w||_2"]
Label(fig[2, 1, Top()], L" \times 10^{-18} ", fontsize = 20.0, halign = :left)
Label(fig[2, 2, Top()], L" \times 10^{-8} ",  fontsize = 20.0, halign = :left)
 title = L"(\varrho, \varrho v, \varrho E)"
#Legend(fig[1, 1:2], handles_colors, labels_colors, title, titlesize = 30, orientation = :horizontal, labelsize = 25.0)
Legend(fig[1, 1:2], ax1, title, titlesize = 30, orientation = :horizontal, labelsize = 25.0)

save(pwd() * "/test_cases/balance/plots/balance_rhoE.pdf", fig)

dd = load_processed_data("well_balancing_Double64_rhotheta.dat", folder_path)
ff = load_processed_data("well_balancing_Float64_rhotheta.dat", folder_path)
ff_theta = load_processed_data("well_balancing_theta_Float64_rhotheta.dat", folder_path)
dd_theta = load_processed_data("well_balancing_theta_Double64_rhotheta.dat", folder_path)
colors = Makie.wong_colors()
stepsize = 2500
stepmarker = 40000
fig = Figure(size = (1200, 600))
kwargs = (xlabel = L"t", xlabelsize = 24, xticklabelsize = 17.0, yticklabelsize = 17.0,  titlesize = 20)
ax1 = Axis(fig[2, 1]; kwargs..., ylabelsize = 25, title = L"\text{Double64}", limits = ((0 ,5000) ,(-0.25, 3) ))
n = 500000
markers = (:circle,)
lines!(ax1, dd[2:n,1], dd[2:n,end-1]./1e-16; label = L"$||u||_2 (T = \text{const})$", linewidth = 2.0, color = colors[1])
lines!(ax1, dd[2:n,1], dd[2:n,end]./1e-16; label = L"$||w||_2 (T = \text{const})$", linewidth = 4.0, color = colors[1], linestyle =:dash)
lines!(ax1, dd_theta[2:n,2], dd_theta[2:n,end-1]./1e-16; label = L"$||u||_2 (\theta = \text{const})$", linewidth = 2.0, color = colors[2])
lines!(ax1, dd_theta[2:stepsize:n,2], dd_theta[2:stepsize:n,end]./1e-16; label = L"$||w||_2 (\theta = \text{const})$", linewidth = 4.0, linestyle =:dash, color = colors[2],)

ax2 = Axis(fig[2, 2]; kwargs..., ylabelsize = 25, title = L"\text{Float64}", limits = ((0 ,5000) ,(-0.25, 6.5) ))
lines!(ax2, ff[2:n,2], ff[2:n,end-1]./1e-9; label = L"$||u||_2$", linewidth = 2.0, color = colors[1])
lines!(ax2, ff[2:n,2], ff[2:n,end]./1e-9; label = L"$||w||_2$", linewidth = 4.0, color = colors[1], linestyle =:dash)
lines!(ax2, ff_theta[2:n,2], ff_theta[2:n,end-1]./1e-9; label = L"$||u||_2$", linewidth = 2.0, color = colors[2])
lines!(ax2, ff_theta[2:n,2], ff_theta[2:n,end]./1e-9; label = L"$||w||_2$", linewidth = 4.0, color = colors[2], linestyle =:dash)

handles_colors = [
    LineElement(color = colors[1], linewidth = 3),
    LineElement(color = colors[2], linewidth = 3), LineElement(color = :black, linewidth = 3, linestyle = :solid),
    LineElement(color = :black, linewidth = 3, linestyle = :dash)
]
labels_colors = [L"T = \text{const}", L"\theta = \text{const}", L"||u||_2", L"||w||_2"]
Label(fig[2, 1, Top()], L" \times 10^{-16} ", fontsize = 20.0, halign = :left)
Label(fig[2, 2, Top()], L" \times 10^{-9} ",  fontsize = 20.0, halign = :left)
 title = L"(\varrho, \varrho v, \varrho \theta)"

#Legend(fig[1, 1:2], handles_colors, labels_colors, title, titlesize = 30, orientation = :horizontal, labelsize = 25.0)
Legend(fig[1, 1:2], ax1, title, titlesize = 30, orientation = :horizontal, labelsize = 25.0)

save(pwd() * "/test_cases/balance/plots/balance_rhotheta.pdf", fig)
