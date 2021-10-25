function fig_print(x_data, y_data, fig_filename = [], 
					legend_string = nothing, xlabel_string = nothing, ylabel_string = nothing, fig_handle = nothing;  
					lw_val = 1.5, N_markers = 10, mks_val = 4, 
					gfs_val = 9, lfs_val = 8, ar_val = :auto, save_file = 1)
	"""
	Figure generation and exporting to file

	fig_print(x_data, y_data, fig_filename, legend_string, xlabel_string, ylabel_string)

	To just save an existing plot "f" to file without modification, use the following code:

	fig_print([], [], fig_filename; fig_handle = f)

	-------------------
	x_data: 1D/2D array, columns are series
	y_data: 2D array, 	 columns are series
	legend_string: row vector (1D array)

	if     fig_handle is given  , then fig_print updates an existing plot
	elseif fig_handle is nothing, then fig_print generates a new plot  

	if 	   x_data or y_data is empty, then fig_print still controls axes, legends, file saving
	"""

	#------ linespec
	linestyle_preset   = repeat([:solid	:solid	:solid	:solid	:dashdot :dash	:dot], 1, 5)
	# choose from [:auto, :solid, :dash, :dot, :dashdot, :dashdotdot]
	markershape_preset = repeat([:circle	:rect	:cross	:none	], 1, 5)
	# choose from [:none, :auto, :circle, :rect, :star5, :diamond, :hexagon, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :heptagon, :octagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x].

	if ~isempty(x_data) && ~isempty(y_data)
	#------ marker index
		marker_idx = zeros(size(x_data, 1))
		marker_idx[round.(Int, LinRange(1, size(x_data, 1), N_markers))] = mks_val*ones(N_markers)
		# marker_idx = zeros(length(fig_handle.series_list[1].plotattributes[:x]))
		# marker_idx[round.(Int, LinRange(1,length(fig_handle.series_list[1].plotattributes[:x]), N_markers))] = mks_val*ones(N_markers)
			
	#------ figure generation / augmentation
		if isnothing(fig_handle)
		fig_handle = plot(x_data, y_data, linewidth = lw_val, linestyle = linestyle_preset,
					label = legend_string, legend = :best,
					markersize = marker_idx, markershape = markershape_preset)
		else
			plot!(fig_handle, x_data, y_data, linewidth = lw_val, linestyle = linestyle_preset[:,length(fig_handle.series_list)+1:end],
					label = legend_string, legend = :best,
					markersize = marker_idx, markershape = markershape_preset[:,length(fig_handle.series_list)+1:end])
		end
	end


	if ~isnothing(fig_handle)
	#------ axes
		plot!(fig_handle, aspect_ratio = ar_val)
		plot!(fig_handle, guidefontsize = gfs_val)
		if ~isnothing(xlabel_string) && ~isempty(xlabel_string)
			plot!(fig_handle, xlabel = xlabel_string)
		end
		if ~isnothing(ylabel_string) && ~isempty(ylabel_string)
			plot!(fig_handle, ylabel = ylabel_string)
		end
			
	#------ legend
		if ~isnothing(legend_string) && length(legend_string) > 1 && typeof(legend_string) != String
			for i in 1:length(fig_handle.series_list)
				fig_handle.series_list[i].plotattributes[:label] = legend_string[i]
			end
		end
		plot!(fig_handle, legendfontsize = lfs_val)

	#------ save
		if save_file == 1
			fig_dir = "Figures"
			mkpath(fig_dir)
			savefig(fig_handle, joinpath(fig_dir, string("Fig_", fig_filename, ".pdf")))
		end
	end

	return fig_handle
end


function equal_AR_3D(x_data, y_data, z_data, legend_string = :false, fig_handle = nothing)
	x12, y12, z12 = extrema(x_data), extrema(y_data), extrema(z_data)
	d = maximum([diff([x12...]),diff([y12...]),diff([z12...])])[1] / 2
	(xm, ym, zm) = (x12, y12, z12) .|> x->(x[1]+x[2])/2

	lw_val, N_markers, mks_val = 1.5, 10, 4
	#------ linespec
	linestyle_preset   = repeat([:solid	:solid	:solid	:solid	:dashdot :dash	:dot], 1, 5)
	markershape_preset = repeat([:circle	:rect	:cross	:none	:none	 :none	:none], 1, 5)
	#------ marker index
	marker_idx = zeros(size(x_data, 1))
	marker_idx[round.(Int, LinRange(1, size(x_data, 1), N_markers))] = mks_val*ones(N_markers)		

	#------ figure generation / augmentation
	if ~isempty(x_data) && ~isempty(y_data) && ~isempty(z_data)
		if isnothing(fig_handle)			
			fig_handle = plot(x_data, y_data, z_data, linewidth = lw_val, linestyle = linestyle_preset,
							label = legend_string, legend = :best,
							markersize = marker_idx, markershape = markershape_preset)
		else
			plot!(fig_handle, x_data, y_data, z_data, linewidth = lw_val, linestyle = linestyle_preset[:,length(fig_handle.series_list)+1:end],
					label = legend_string, legend = :best,
					markersize = marker_idx, markershape = markershape_preset[:,length(fig_handle.series_list)+1:end])

		end
	end
	plot!(fig_handle, aspect_ratio = :equal, xlims=(xm-d,xm+d), ylims=(ym-d,ym+d), zlims=(zm-d,zm+d))
	
	return fig_handle
end