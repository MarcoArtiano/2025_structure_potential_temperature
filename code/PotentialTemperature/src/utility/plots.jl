using LaTeXStrings

# function load_processed_data(file_name, folder_path)
# 	file_path = joinpath(folder_path, file_name)
# 	return readdlm(file_path)
# end

function create_plot!(data_tuple, label_tuple, ylabel, index; position = :lt)

	fig = Figure()
	ax = Axis(fig[1, 1]; xlabel = L"t", ylabel = ylabel)
	#xlims!(0, 40)
	#ylims!(-1e-5, 1e-5)

	for (var, label) in zip(data_tuple, label_tuple)
		plot_var!(var, index, label, ax)
	end
	axislegend(position = position)
	return fig
end

function plot_var!(var, index, legend, ax)

	lines!(ax, var[:,1], var[:, index]; label = legend)

end

function create_plot2!(data_tuple, label_tuple, ax, index)

	for (var, label) in zip(data_tuple, label_tuple)
		plot_var!(var, index, label, ax)
	end
end