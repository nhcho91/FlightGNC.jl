function fig_print(x_data, y_data, fig_filename = [], 
					legend_string = nothing, xlabel_string = nothing, ylabel_string = nothing, fig_handle = nothing;  
					lw_val = 1.5, lgnd_val = :best, N_markers = 10, mks_val = 4, 
					gfs_val = 12, lfs_val = 8, ar_val = :auto, save_file = 1, kwargs...)
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
	#------ figure generation / augmentation
		if isnothing(fig_handle)
		fig_handle = plot(x_data, y_data, linewidth = lw_val, linestyle = linestyle_preset,
					label = legend_string, legend = lgnd_val)
					# markersize = marker_idx, markershape = markershape_preset)

			#--- marker index (post-processing)
			for i in 1 : length(fig_handle.series_list)
				marker_idx = zeros(length(fig_handle.series_list[i].plotattributes[:x]))
				marker_idx[round.(Int, LinRange(1,length(fig_handle.series_list[i].plotattributes[:x]), N_markers))] = mks_val*ones(N_markers)
				fig_handle.series_list[i].plotattributes[:markersize]  = marker_idx
				fig_handle.series_list[i].plotattributes[:markershape] = markershape_preset[i]
			end

		else
			#--- marker index (pre-processing)
			marker_idx = zeros(size(x_data, 1))
			marker_idx[round.(Int, LinRange(1, size(x_data, 1), N_markers))] = mks_val*ones(N_markers)
			plot!(fig_handle, x_data, y_data, linewidth = lw_val, linestyle = linestyle_preset[:,length(fig_handle.series_list)+1:end],
					label = legend_string, legend = lgnd_val,
					markersize = marker_idx, markershape = markershape_preset[:,length(fig_handle.series_list)+1:end])
		end
	end


	if ~isnothing(fig_handle)
		plot!(fig_handle, kwargs...)
		
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

		# Legend location setting should precede legend fontsize setting
		plot!(fig_handle, legend = lgnd_val)
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


# view_result() is developed for use with CTPG_train()
function view_result(i_plot_list, fwd_ensemble_sol, loss_history, save_fig = true;
    vars_x = nothing, x_names = [], vars_u = nothing, u_names = [], vars_y = nothing, y_names = [], vars_y_NN = nothing, y_NN_names = [], plot_kwargs...)

    # labels, variables to be plotted, and layouts setting
    (dim_x_aug, dim_u, dim_y, dim_y_NN) = (length(fwd_ensemble_sol[1].sol[1]), length(fwd_ensemble_sol[1].u[1]), length(fwd_ensemble_sol[1].y[1]), length(fwd_ensemble_sol[1].y_NN[1]))
    (xlabel_t, ylabel_x, ylabel_u, ylabel_y, ylabel_y_NN) = ("\$t\$", "state", "input", "output", "NN output")

    if ~isnothing(vars_x)
        layout_x = (length(vars_x), 1)
        if ~isempty(x_names)
            ylabel_x = x_names[vars_x']
        end
    else
        layout_x = (dim_x_aug - 1, 1)
        vars_x = 1:dim_x_aug-1
    end

    if ~isnothing(vars_u)
        layout_u = (length(vars_u), 1)
        if ~isempty(u_names)
            ylabel_u = u_names[vars_u']
        end
    else
        layout_u = (dim_u, 1)
        vars_u = 1:dim_u
    end

    if ~isnothing(vars_y)
        layout_y = (length(vars_y), 1)
        if ~isempty(y_names)
            ylabel_y = y_names[vars_y']
        end
    else
        layout_y = (dim_y, 1)
        vars_y = 1:dim_y
    end

    if ~isnothing(vars_y_NN)
        layout_y_NN = (length(vars_y_NN), 1)
        if ~isempty(y_NN_names)
            ylabel_y_NN = y_NN_names[vars_y_NN']
        end
    else
        layout_y_NN = (dim_y_NN, 1)
        vars_y_NN = 1:dim_y_NN
    end

    # plotting
    if isempty(i_plot_list)
        i_plot_list = 1:length(fwd_ensemble_sol)
    end

    (f_x, f_u, f_y, f_y_NN) = (plot(), plot(), plot(), plot())
    for i in i_plot_list
        # data preprocessing
        @unpack sol, u, y, y_NN = fwd_ensemble_sol[i]
        @unpack t = sol

        u = reduce(vcat, u')     # hcat(u...)'
        y = reduce(vcat, y')     # hcat(y...)'
        y_NN = reduce(vcat, y_NN')  # hcat(y_NN...)'

        if i == first(i_plot_list)
            f_x = plot(sol, vars = vars_x, layout = layout_x, label = :false, xlabel = xlabel_t, ylabel = ylabel_x, size = (800, 160 * length(vars_x)); plot_kwargs...)

            f_u = plot(t, u[:, vars_u], layout = layout_u, label = :false, xlabel = xlabel_t, ylabel = ylabel_u; plot_kwargs...)

            f_y = plot(t, y[:, vars_y], layout = layout_y, label = :false, xlabel = xlabel_t, ylabel = ylabel_y; plot_kwargs...)

            f_y_NN = plot(t, y_NN[:, vars_y_NN], layout = layout_y_NN, label = :false, xlabel = xlabel_t, ylabel = ylabel_y_NN; plot_kwargs...)
        else
            plot!(f_x, sol, vars = vars_x, layout = layout_x, label = :false, xlabel = xlabel_t, ylabel = ylabel_x, size = (800, 160 * length(vars_x)); plot_kwargs...)

            plot!(f_u, t, u[:, vars_u], layout = layout_u, label = :false, xlabel = xlabel_t, ylabel = ylabel_u; plot_kwargs...)

            plot!(f_y, t, y[:, vars_y], layout = layout_y, label = :false, xlabel = xlabel_t, ylabel = ylabel_y; plot_kwargs...)

            plot!(f_y_NN, t, y_NN[:, vars_y_NN], layout = layout_y_NN, label = :false, xlabel = xlabel_t, ylabel = ylabel_y_NN; plot_kwargs...)
        end
    end
    f_L = plot(loss_history, label = :false, xlabel = "iteration", ylabel = "\$L\$"; plot_kwargs...)

    display(f_x)
    display(f_u)
    display(f_y)
    display(f_y_NN)
    display(f_L)

    if save_fig
        savefig(f_x,    "f_x.pdf")
        savefig(f_u,    "f_u.pdf")
        savefig(f_y,    "f_y.pdf")
        savefig(f_y_NN, "f_y_NN.pdf")
        savefig(f_L,    "f_L.pdf")
    end

    return f_x, f_u, f_y, f_y_NN, f_L
end
