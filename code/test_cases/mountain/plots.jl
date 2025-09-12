  using JLD2
  using PotentialTemperature
  using PotentialTemperature.Trixi
  using CairoMakie
  @load "mountain.jld2" sol semi cells_per_dimension equations
  polydeg = 5
    x, y, data_eulernc = ContourData(sol.u[end], semi, cells_per_dimension, equations)
    # Per u tra -2 e 2 con passo 0.2
    up = data_eulernc[2, :, :] ./ data_eulernc[1, :, :] .-10.0
    
    wp = data_eulernc[3, :, :] ./ data_eulernc[1, :, :]
# Per w tra -0.5 e 0.5 con passo 0.05
    maxu = maximum(up)
    minu = minimum(up)
    maxw = maximum(wp)
    minw = minimum(wp)
    levels_u = collect(minu:0.2:maxu)
    levels_w = collect(minw:0.05:maxw)
    h = Figure(size = (850, 400))
    labelsize = 20
    levels = 13
    kwargs = (xlabel = L"$x$ [km]", xlabelsize = labelsize, ylabelsize = labelsize, limits = ((-10, 10), (0, 10)), xticklabelsize = 17.0, yticklabelsize = 17.0)

    Axis(h[1, 1]; kwargs..., ylabel = L"$z$ [km]")

    c = contour!(x ./ 1e3, y ./ 1e3, up, levels = levels_u, color = :black)
    Axis(h[1, 2]; kwargs...)

    c = contour!(x ./ 1e3, y ./ 1e3, wp , levels = levels_w, color = :black)
    save(pwd() * "/test_cases/mountain/plots/contour_mountain_new_$(cells_per_dimension[1])x$(cells_per_dimension[2])_CFL1_polydeg$(polydeg)_u_w.pdf", h)
