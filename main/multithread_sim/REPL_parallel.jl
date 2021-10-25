##
include("main_sim_parallel.jl")


"""
# ----------------------------------------------------------------------------------------------
## simulation with fixed parameters (single run)
# χ_M_0 = deg2rad(60)
# v_M_0 = V_M_0 * [sin(χ_M_0); cos(χ_M_0)]

# χ_f_d = deg2rad(180)
# s_BPNG.v̂_f_d = [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)]
# df, f = sim_plot(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG, A_M_max);
# display(f)


# ----------------------------------------------------------------------------------------------
## simulation with one parameter variation
# Variable 1
χ_M_0_list  = deg2rad.(30:5:90)
v_M_0_list  = χ_M_0_list |> Map(χ_M_0 -> V_M_0 * [sin(χ_M_0); cos(χ_M_0)]) |> collect

# Variable 2
χ_f_d_list  = deg2rad.(90:5:270)
v̂_f_d_list  = χ_f_d_list    |> Map(χ_f_d -> [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)])             |> collect
s_BPNG_list = v̂_f_d_list    |> Map(v̂_f_d -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget_2D, s_Bias))    |> collect

# size(df_list)

# df_list     = v_M_0_list |> Map(v_M_0 -> main( p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max = A_M_max) ) |> tcollect
df_list     = s_BPNG_list   |> Map(s_BPNG -> main( p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max = A_M_max))          |> tcollect

t_list      = df_list   |> Map(df -> df.time)     |> collect
# sol_list    = df_list   |> Map(df -> df.sol)      |> collect
x_list      = df_list   |> Map(df -> df.sol |> Map(datum -> datum.pursuer.p[1])     |> collect)    |> collect
y_list      = df_list   |> Map(df -> df.sol |> Map(datum -> datum.pursuer.p[2])     |> collect)    |> collect
χ_list      = df_list   |> Map(df -> df.sol |> Map(datum -> atan(datum.pursuer.v[1], datum.pursuer.v[2]) )     |> collect)    |> collect
A_list      = df_list   |> Map(df -> df.sol |> Map(datum -> norm(datum.pursuer.u))  |> collect)    |> collect
r_list      = df_list   |> Map(df -> df.sol |> Map(datum -> datum.r)                |> collect)    |> collect
σ_M_list    = df_list   |> Map(df -> df.sol |> Map(datum -> datum.σ_M)              |> collect)    |> collect
miss_list   = r_list    |> Map(r  -> r[end])    |> collect

f_xy = plot(size = (800,600))
for i in 1:length(t_list)
    plot!(f_xy, x_list[i], y_list[i], label = :false, xlabel = "x [m]", ylabel = "y [m]")
end
f_xy = fig_print([], [], "xy"; fig_handle = f_xy, ar_val = :equal) 	# To just save the result
display(f_xy)

f_r = plot()
for i in 1:length(t_list)
    plot!(f_r, t_list[i], r_list[i], label = :false, xlabel = "t [s]", ylabel = "r [m]")
end
f_r = fig_print([], [], "r"; fig_handle = f_r) 	# To just save the result
display(f_r)

# f_miss = fig_print(rad2deg.(χ_f_d_list), miss_list, "miss", :false, "χ_f_d [deg]", "r_f [m]")
# display(f_miss)

f_σ = plot()
for i in 1:length(t_list)
    plot!(f_σ, t_list[i], rad2deg.(σ_M_list[i]), label = :false, xlabel = "t [s]", ylabel = "σ_M [deg]")
end
f_σ = fig_print([], [], "sigma_M"; fig_handle = f_σ) 	# To just save the result
display(f_σ)

f_χ = plot()
for i in 1:length(t_list)
    plot!(f_χ, t_list[i], rad2deg.(χ_list[i]), label = :false, xlabel = "t [s]", ylabel = "χ_M [deg]")
    plot!(f_χ, t_list[i][[1,end],1], rad2deg.(χ_f_d_list[i]*ones(2)), label = :false, linestyle = :dot, linewidth = 0.5)
end
f_χ = fig_print([], [], "chi_M"; fig_handle = f_χ) 	# To just save the result
display(f_χ)

f_A = plot()
for i in 1:length(t_list)
    plot!(f_A, t_list[i], A_list[i], label = :false, xlabel = "t [s]", ylabel = "A_M [m/s^2]")
end
f_A = fig_print([], [], "A_M"; fig_handle = f_A) 	# To just save the result
display(f_A)



# ----------------------------------------------------------------------------------------------
## simulation with two parameters variation
df_list_list  = v_M_0_list |> Map( v_M_0 -> s_BPNG_list |> Map(s_BPNG -> main( p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max = A_M_max)) |> tcollect  ) |> tcollect

t_list      = df_list_list |> Map( df_ic -> df_ic   |> Map(df_param -> df_param.time)      |> collect) |> collect
# sol_list    = df_list_list |> Map( df_ic -> df_ic   |> Map(df_param -> df_param.sol)       |> collect) |> collect

x_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> datum.pursuer.p[1])                              |> collect  )  |> collect )    |> collect
y_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> datum.pursuer.p[2])                              |> collect  )  |> collect )    |> collect
χ_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> atan(datum.pursuer.v[1], datum.pursuer.v[2]) )   |> collect  )  |> collect )    |> collect
A_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> norm(datum.pursuer.u))                           |> collect  )  |> collect )    |> collect
r_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> datum.r)                                         |> collect  )  |> collect )    |> collect
σ_M_list    = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> datum.σ_M)                                       |> collect  )  |> collect )    |> collect
miss_list   = r_list       |> Map( r_ic_list -> r_ic_list   |> Map(r  -> r[end])    |> collect ) |> collect

miss_list   = hcat(miss_list...)'

#
f_miss_surf = contour(rad2deg.(χ_M_0_list), rad2deg.(χ_f_d_list), miss_list', xlabel = "χ_M_0 [deg]", ylabel = "χ_f_d [deg]", fill = true)
# plot(rad2deg.(χ_M_0_list), rad2deg.(χ_f_d_list), miss_list, st = [:surface, :contourf], layout = 2)
f_miss_surf = fig_print([], [], "miss"; fig_handle = f_miss_surf) 	# To just save the result

i = 2 # in 1:length(v_M_0_list)

for i in 1:length(v_M_0_list)
    f_xy = plot(size = (800,600))

    for j in 1:length(s_BPNG_list)
        plot!(f_xy, x_list[i][j], y_list[i][j], label = :false, xlabel = "x [m]", ylabel = "y [m]")
    end
    f_xy = fig_print([], [], "xy_ic_"*string(i); fig_handle = f_xy, ar_val = :equal) 	# To just save the result
    display(f_xy)
end

# f_r = plot()
# for j in 1:length(s_BPNG_list)
#     plot!(f_r, t_list[i][j], r_list[i][j], label = :false, xlabel = "t [s]", ylabel = "r [m]")
# end
# f_r = fig_print([], [], "r_ic_"*string(i); fig_handle = f_r) 	# To just save the result
# display(f_r)

# f_σ = plot()
# for j in 1:length(s_BPNG_list)
#     plot!(f_σ, t_list[i][j], rad2deg.(σ_M_list[i][j]), label = :false, xlabel = "t [s]", ylabel = "σ_M [deg]")
# end
# f_σ = fig_print([], [], "sigma_M_ic_"*string(i); fig_handle = f_σ) 	# To just save the result
# display(f_σ)

# f_χ = plot()
# for j in 1:length(s_BPNG_list)
#     plot!(f_χ, t_list[i][j], rad2deg.(χ_list[i][j]), label = :false, xlabel = "t [s]", ylabel = "χ_M [deg]")
#     plot!(f_χ, t_list[i][j][[1,end],1], rad2deg.(χ_f_d_list[j]*ones(2)), label = :false, linestyle = :dot, linewidth = 0.5)
# end
# f_χ = fig_print([], [], "chi_M_ic_"*string(i); fig_handle = f_χ) 	# To just save the result
# display(f_χ)

# f_A = plot()
# for j in 1:length(s_BPNG_list)
#     plot!(f_A, t_list[i][j], A_list[i][j], label = :false, xlabel = "t [s]", ylabel = "A_M [m/s^2]")
# end
# f_A = fig_print([], [], "A_M_ic_"*string(i); fig_handle = f_A) 	# To just save the result
# display(f_A)



# ---------------------------------------------------------------------------------------------------
## 3D guidance
include("main_missiles.jl")
include("sim_missiles.jl")


##
# Parallel execution (works only in REPL mode)
# Threads.nthreads()

dim = 3

γ_f_d_list  = deg2rad.(-90:5:0)
χ_f_d_list  = deg2rad.(90:5:270)
v̂_f_d_list  = γ_f_d_list    |> Map( γ_f_d -> χ_f_d_list |> Map(χ_f_d -> [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)])             |> collect ) |> collect
s_BPNG_list = v̂_f_d_list    |> Map( v̂_f_d_γ -> v̂_f_d_γ  |> Map(v̂_f_d -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))       |> collect ) |> collect

df_list_list  = s_BPNG_list |> Map( s_BPNG_γ -> s_BPNG_γ |> Map(s_BPNG -> main( p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max = A_M_max)) |> tcollect  ) |> tcollect

x_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> datum.pursuer.p[1])  |> collect  )  |> collect )    |> collect
y_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> datum.pursuer.p[2])  |> collect  )  |> collect )    |> collect
z_list      = df_list_list |> Map( df_ic -> df_ic  |> Map(df_param -> df_param.sol |> Map(datum -> datum.pursuer.p[3])  |> collect  )  |> collect )    |> collect

f_3D = plot(aspect_ratio = :equal)
for i in 1:length(γ_f_d_list)
    for j in 1:length(χ_f_d_list)
        plot!(f_3D, x_list[i][j], y_list[i][j], z_list[i][j], label=:false)
    end
end
plot!(f_3D, xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]", camera = (45,30), size = (800,800), 
xlims = (0,5.5E3), ylims = (-1E3, 4.5E3), zlims = (-1E3,4E3))
f_3D = fig_print([], [], "Trajectory_3D"; fig_handle = f_3D, ar_val = :equal) 	# To just save the result

"""


## ---------------------------------------------------------------------------------------------------
# massive numerical analysis example
# ---------------------------------------------------------------------------------------------------
using LaTeXStrings
include("main_sim_parallel.jl")


# ---------------------------------------------------------------------------------------------------
## List of parameters
# Variable 1
γ_xref_M_0_list = deg2rad.(0:5:60)
χ_M_0_list  = pi / 2 .- γ_xref_M_0_list # deg2rad.(30:5:90)
v_M_0_list  = χ_M_0_list |> Map(χ_M_0 -> V_M_0 * [sin(χ_M_0); cos(χ_M_0)]) |> collect

# Variable 2
γ_xref_f_d_list = deg2rad.(-180:5:0)
χ_f_d_list  = pi / 2 .- γ_xref_f_d_list # deg2rad.(90:5:270)
v̂_f_d_list  = χ_f_d_list    |> Map(χ_f_d -> [cos(γ_f_d) * sin(χ_f_d); cos(γ_f_d) * cos(χ_f_d); sin(γ_f_d)])             |> collect
# s_BPNG_list = v̂_f_d_list    |> Map(v̂_f_d -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget_2D, s_Bias))    |> collect

# Variable 3-4
α_list      = 0.6:0.2:1
n_list      = 0:1:3
# s_BPNG_list   = α_list      |> Map( α -> n_list  |> Map( n -> begin 
#                     s_Bias  = ComponentArray(α = α, δ = 0.01, n = n, r_ref = 10E3, k = 4, m = 10);
#                     BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias)
#                 end) |> collect ) |> collect

main(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max=A_M_max)  
# --------> the above line is to execute main() once for compilation before passing to multiple threads

# ---------------------------------------------------------------------------------------------------
## simulation for paper: gain design parameter variation
mutable struct Results
    α
    n
    r_f_result
    r_min_result
end
s_sim = Array{Results}(undef, length(α_list), length(n_list))

for i_α in 1:length(α_list)
    for i_n in 1:length(n_list)
        α       = α_list[i_α]
        n       = n_list[i_n]
        s_Bias  = ComponentArray(α=α, δ=0.01, n=n, r_ref=5E3, k=4, m=10)
        s_BPNG_list = v̂_f_d_list    |> Map(v̂_f_d -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget_2D, s_Bias))    |> collect


        # execution
        df_list_list  = v_M_0_list |> Map(v_M_0 -> s_BPNG_list |> Map(s_BPNG -> main(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max=A_M_max)) |> tcollect) |> tcollect
        
        # data collection
        t_list      = df_list_list |> Map(df_ic -> df_ic   |> Map(df_param -> df_param.time)      |> collect) |> collect
        # sol_list    = df_list_list |> Map( df_ic -> df_ic   |> Map(df_param -> df_param.sol)       |> collect) |> collect
        x_list      = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> datum.pursuer.p[1])                              |> collect)  |> collect)    |> collect
        y_list      = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> datum.pursuer.p[2])                              |> collect)  |> collect)    |> collect
        χ_list      = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> atan(datum.pursuer.v[1], datum.pursuer.v[2]))    |> collect)  |> collect)    |> collect
        γ_xref_list = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> atan(datum.pursuer.v[2], datum.pursuer.v[1]))    |> collect)  |> collect)    |> collect
        λ_xref_list = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> atan(datum.evador.p[2] - datum.pursuer.p[2], datum.evador.p[1] - datum.pursuer.p[1]))    |> collect)  |> collect)    |> collect
        A_list      = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> norm(datum.pursuer.u))                           |> collect)  |> collect)    |> collect
        r_list      = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> datum.r)                                         |> collect)  |> collect)    |> collect
        σ_M_list    = df_list_list |> Map(df_ic -> df_ic  |> Map(df_param -> df_param.sol  |> Map(datum -> datum.σ_M)                                       |> collect)  |> collect)    |> collect
        
        r_f_result   = r_list       |> Map(r_ic_list -> r_ic_list   |> Map(r  -> r[end])    |> collect) |> collect
        r_min_result = r_list       |> Map(r_ic_list -> r_ic_list   |> Map(r  -> minimum(r)) |> collect) |> collect
       
        σ_M_signed_list  = γ_xref_list .- λ_xref_list
        γ_xref_pred_list = γ_xref_list .- N / (N - 1) .* σ_M_signed_list

        r_f_result      = hcat(r_f_result...)'
        r_min_result    = hcat(r_min_result...)'
        s_sim[i_α, i_n] = Results(α, n, r_f_result, r_min_result)

        γ_xref_M_0_list = 90 .- rad2deg.(χ_M_0_list)
        γ_xref_f_d_list = 90 .- rad2deg.(χ_f_d_list)
        K_r_max         = (N - 1 + 0.01) * 2^n
        γ_xref_f_d_ub = min(rad2deg.(((N - 1) * tan(σ_M_lim) / K_r_max)^(1 / α)) .- γ_xref_M_0_list / (N - 1), zeros(length(γ_xref_M_0_list)))
        γ_xref_f_d_lb =      rad2deg.(-((N - 1) * tan(σ_M_lim) / K_r_max)^(1 / α)) .- γ_xref_M_0_list / (N - 1)
        
        # r_f
        f_r_f = contour(γ_xref_M_0_list, γ_xref_f_d_list, r_f_result', 
                        xlabel=L"\gamma_0 \textrm{~[deg]}", ylabel=L"\gamma_f_d \textrm{~[deg]}", fill=true, seriescolor=:coolwarm, ylims=extrema(γ_xref_f_d_list))
        plot!(f_r_f, γ_xref_M_0_list, [γ_xref_f_d_lb γ_xref_f_d_ub], linecolor=:red, linestyle=:dash, linewidth=2,
                    label=["LB" "UB"], ylims=extrema(γ_xref_f_d_list)) # title="r_f: α = $α, n = $n"
        f_r_f = fig_print([], [], "r_f_alpha_$(i_α)_n_$(i_n)", [], [], [], f_r_f) 	# To just save the result
        display(f_r_f)
        
        # r_min
        f_r_min = contour(γ_xref_M_0_list, γ_xref_f_d_list, r_min_result', 
        xlabel=L"\gamma_0 \textrm{~[deg]}", ylabel=L"\gamma_f_d \textrm{~[deg]}", fill=true, seriescolor=:coolwarm, ylims=extrema(γ_xref_f_d_list))
        plot!(f_r_min, γ_xref_M_0_list, [γ_xref_f_d_lb γ_xref_f_d_ub], linecolor=:red, linestyle=:dash, linewidth=2,
                    label=["LB" "UB"], ylims=extrema(γ_xref_f_d_list)) # title="r_f: α = $α, n = $n"
        f_r_min = fig_print([], [], "r_min_alpha_$(i_α)_n_$(i_n)", [], [], [], f_r_min) 	# To just save the result
        display(f_r_min)

        # OX
        f_OX = plot()
        for i in 1:length(γ_xref_M_0_list)
            for j in 1:length(γ_xref_f_d_list)
                if r_min_result[i,j] < 4 && r_f_result[i,j] < 4
                    plot!(f_OX, γ_xref_M_0_list[i] * ones(2), γ_xref_f_d_list[j] * ones(2), markershape=:circle, markercolor=:blue, label=:false)
                else
                    plot!(f_OX, γ_xref_M_0_list[i] * ones(2), γ_xref_f_d_list[j] * ones(2), markershape=:xcross, markercolor=:red, label=:false)
                end
            end
        end
        plot!(f_OX, xlabel=L"\gamma_0 \textrm{~[deg]}", ylabel=L"\gamma_f_d \textrm{~[deg]}")
        f_OX = fig_print([], [], "OX_alpha_$(i_α)_n_$(i_n)", [], [], [], f_OX) 	# To just save the result
        display(f_OX)

        # e_γ_f
        for i in 1:length(γ_xref_M_0_list)
            f_e_γ_f = plot!()

            for j in 1:length(γ_xref_f_d_list)    
                γ_xref_f_d = γ_xref_f_d_list[j]
                e_γ_f_list = 180 / pi .* γ_xref_pred_list[i][j] .- γ_xref_f_d
                plot!(f_e_γ_f, t_list[i][j], abs.(e_γ_f_list), label=:false, xlabel="t [s]", ylabel="|e_γ_f| [deg]", ylims=(0, 60))
            end
            f_e_γ_f = fig_print([], [], "e_gamma_f_alpha_$(i_α)_n_$(i_n)_ic_$i", [], [], [], f_e_γ_f) 	# To just save the result
            # display(f_e_γ_f)
        end
        

        # xy
        for i in 1:length(v_M_0_list)
            f_xy = plot(size=(800, 600))

            for j in 1:length(s_BPNG_list)
                plot!(f_xy, x_list[i][j], y_list[i][j], label=:false, xlabel="x [m]", ylabel="y [m]", ylims=(-500, 4.5E3))
            end
            f_xy = fig_print([], [], "xy_alpha_$(i_α)_n_$(i_n)_ic_$i", [], [], [], f_xy, ar_val=:equal) 	# To just save the result
            # display(f_xy)
        end

        # f_r = plot()
        # for j in 1:length(s_BPNG_list)
        #     plot!(f_r, t_list[i][j], r_list[i][j], label = :false, xlabel = "t [s]", ylabel = "r [m]")
        # end
        # f_r = fig_print([], [], "r", [], [], [], f_r) 	# To just save the result
        # display(f_r)


        # f_σ = plot()
        # for j in 1:length(s_BPNG_list)
        #     plot!(f_σ, t_list[i][j], rad2deg.(σ_M_list[i][j]), label = :false, xlabel = "t [s]", ylabel = "σ_M [deg]")
        # end
        # f_σ = fig_print([], [], "sigma_M", [], [], [], f_σ) 	# To just save the result
        # display(f_σ)

        # f_χ = plot()
        # for j in 1:length(s_BPNG_list)
        #     plot!(f_χ, t_list[i][j], rad2deg.(χ_list[i][j]), label = :false, xlabel = "t [s]", ylabel = "χ_M [deg]")
        #     plot!(f_χ, t_list[i][j][[1,end],1], rad2deg.(χ_f_d_list[j]*ones(2)), label = :false, linestyle = :dot, linewidth = 0.5)
        # end
        # f_χ = fig_print([], [], "chi_M", [], [], [], f_χ) 	# To just save the result
        # display(f_χ)

        # f_A = plot()
        # for j in 1:length(s_BPNG_list)
        #     plot!(f_A, t_list[i][j], A_list[i][j], label = :false, xlabel = "t [s]", ylabel = "A_M [m/s^2]")
        # end
        # f_A = fig_print([], [], "A_M", [], [], [], f_A) 	# To just save the result
        # display(f_A)

    end
end

