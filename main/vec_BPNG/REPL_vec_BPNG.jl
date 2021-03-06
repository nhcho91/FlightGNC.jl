include("main_vec_BPNG.jl")

# Execute main() once for compilation before passing to multiple threads
# df = main(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max=A_M_max)

## [Manual Choice]
case = 6

# List of simulation parameters
if case == 1
    γ_f_d   = deg2rad(-70)
    s_Bias  = ComponentArray(α=1, δ=0, n=1, r_ref=10E3, k=0, m=0, k̂_d=[1; 0; 0], i_Ω_μ = 0, i_σ_M_lim = 0)

    χ_f_d_list  = deg2rad.(0:45:315)
    v̂_f_d_list  = χ_f_d_list    |> Map(χ_f_d -> [cos(γ_f_d) * sin(χ_f_d); cos(γ_f_d) * cos(χ_f_d); sin(γ_f_d)])     |> collect
    s_BPNG_list = v̂_f_d_list    |> Map(v̂_f_d -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))   |> collect

    # γ_f_d_list  = deg2rad.(-90:45:90)
    # χ_f_d_list  = deg2rad.(0:45:180)
    # v̂_f_d_list  = γ_f_d_list    |> Map(γ_f_d -> χ_f_d_list |> Map(χ_f_d -> [cos(γ_f_d) * sin(χ_f_d); cos(γ_f_d) * cos(χ_f_d); sin(γ_f_d)])             |> collect) |> collect
    # s_BPNG_list = v̂_f_d_list    |> Map(v̂_f_d_γ -> v̂_f_d_γ  |> Map(v̂_f_d -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))           |> collect) |> collect

    label_string = Vector{String}(undef, length(s_BPNG_list))
    for i in 1:length(s_BPNG_list)
        label_string[i] = "\$\\chi_{f_{d}} = $(round(Int, rad2deg(χ_f_d_list[i])))\$"
    end
    camera_angle = (50, 40)

elseif case == 2
    γ_f_d   = deg2rad(-70)
    v̂_f_d   = [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)]

    n_list      = 0:1:4
    s_Bias_list = n_list        |> Map(n -> ComponentArray(α=1, δ=0, n=n, r_ref=10E3, k=0, m=0, k̂_d=[1; 0; 0], i_Ω_μ = 0, i_σ_M_lim = 0))      |> collect
    s_BPNG_list = s_Bias_list   |> Map(s_Bias -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))              |> collect

    label_string = Vector{String}(undef, length(s_BPNG_list))
    for i in 1:length(s_BPNG_list)
        label_string[i] = "\$n = $(n_list[i])\$"
    end
    camera_angle = (30, 30)

elseif case == 3
    γ_f_d   = deg2rad(-70)
    χ_f_d   = deg2rad(225)
    v̂_f_d   = [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)]
    s_Bias  = ComponentArray(α=1, δ=0, n=1, r_ref=10E3, k=100, m=0, k̂_d=[1; 0; 0], i_Ω_μ = 0, i_σ_M_lim = 1)

    σ_M_lim_list = deg2rad.(45:5:60)
    s_BPNG_list = σ_M_lim_list   |> Map(σ_M_lim -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))            |> collect

    label_string = Vector{String}(undef, length(s_BPNG_list))
    for i in 1:length(s_BPNG_list)
        label_string[i] = "\$\\sigma_{\\lim} = $(round(Int, rad2deg(σ_M_lim_list[i])))\$"
    end
    camera_angle = (30, 30)

elseif case == 4
    Ω_μ_0_list  = [0; 0.1; 0.15; 0.2; 0.25]
    s_Bias_list = Ω_μ_0_list    |> Map(δ -> ComponentArray(α=1, δ=δ, n=1, r_ref=10E3, k=4, m=3, k̂_d=[1; 0; 0], i_Ω_μ = 2, i_σ_M_lim = 1))      |> collect
    s_BPNG_list = s_Bias_list   |> Map(s_Bias -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))              |> collect

    label_string = Vector{String}(undef, length(s_BPNG_list))
    for i in 1:length(s_BPNG_list)
        label_string[i] = "\$\\Omega_{\\mu_{0}}= $(Ω_μ_0_list[i])\$"
    end
    camera_angle = (30, 45)

elseif case == 5
    k̂_d_list    = [[1; 0; 0], [1 / sqrt(2); 1 / sqrt(2); 0], [0; 1; 0]]
    s_Bias_list = k̂_d_list      |> Map(k̂_d -> ComponentArray(α=1, δ=0, n=1, r_ref=10E3, k=4, m=20, k̂_d=k̂_d, i_Ω_μ = 1, i_σ_M_lim = 1))         |> collect
    s_BPNG_list = s_Bias_list   |> Map(s_Bias -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))              |> collect

    label_string = ["\$\\hat{\\mathbf{k}}_{d} = \\left[1;0;0\\right]\$" 
                    "\$\\hat{\\mathbf{k}}_{d} = \\left[1/\\sqrt{2};1/\\sqrt{2};0\\right]\$" 
                    "\$\\hat{\\mathbf{k}}_{d} = \\left[0;1;0\\right]\$"] 
    camera_angle = (35, 45)

elseif case == 6
    n       = 1.0
    N_σ     = (n + 2.0) * (n + 3.0)
    γ_f_d   = deg2rad(-70)
    s_Bias  = ComponentArray(α=1, δ=0, n=n, r_ref=10E3, k=0, m=3, k̂_d=[1; 0; 0], i_Ω_μ = 0, i_σ_M_lim = 0)

    χ_f_d_list  = deg2rad.(0:45:315)
    v̂_f_d_list  = χ_f_d_list    |> Map(χ_f_d -> [cos(γ_f_d) * sin(χ_f_d); cos(γ_f_d) * cos(χ_f_d); sin(γ_f_d)])     |> collect

    # Quaternion_IACG
    s_BPNG_list = v̂_f_d_list    |> Map(v̂_f_d -> BPNG(N_σ, dim, σ_M_lim, v̂_f_d, Bias_Quaternion_IACG, s_Bias))   |> collect
    # vec_BPNG
    # s_BPNG_list = v̂_f_d_list    |> Map(v̂_f_d -> BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget, s_Bias))   |> collect

    label_string = Vector{String}(undef, length(s_BPNG_list))
    for i in 1:length(s_BPNG_list)
        label_string[i] = "\$\\chi_{f_{d}} = $(round(Int, rad2deg(χ_f_d_list[i])))\$"
    end
    camera_angle = (50, 40)
end


## ---------------------------------------------------------------------------------------------------
# mutable struct Results
#     r_f_result
#     r_min_result
#     t_list
#     x_list
#     y_list
#     A_list
#     r_list
#     σ_M_list    
# end
# s_sim = Array{Results}(undef, length(s_BPNG_list))

# multithread execution
df_list  = s_BPNG_list |> Map(s_BPNG -> main(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max=A_M_max)) |> tcollect
        
# data collection
t_list      = df_list |> Map(df -> df.time)      |> collect
sol_list    = df_list |> Map(df -> df.sol)       |> collect
x_list      = df_list |> Map(df -> df.sol  |> Map(datum -> datum.pursuer.p[1])          |> collect)  |> collect
y_list      = df_list |> Map(df -> df.sol  |> Map(datum -> datum.pursuer.p[2])          |> collect)  |> collect
z_list      = df_list |> Map(df -> df.sol  |> Map(datum -> datum.pursuer.p[3])          |> collect)  |> collect
A_list      = df_list |> Map(df -> df.sol  |> Map(datum -> norm(datum.pursuer.u))       |> collect)  |> collect
r_list      = df_list |> Map(df -> df.sol  |> Map(datum -> datum.r)                     |> collect)  |> collect
σ_M_list    = df_list |> Map(df -> df.sol  |> Map(datum -> datum.σ_M)                   |> collect)  |> collect
e_v̂_f_list  = df_list |> Map(df -> df.sol  |> Map(datum -> datum.e_v̂_f)                 |> collect)  |> collect
A_bias_list = df_list |> Map(df -> df.sol  |> Map(datum -> norm(datum.a_M_bias))        |> collect)  |> collect
Ω_bias_list = df_list |> Map(df -> df.sol  |> Map(datum -> norm(datum.ω_bias))          |> collect)  |> collect
μ_list      = df_list |> Map(df -> df.sol  |> Map(datum -> datum.μ)                     |> collect)  |> collect
Ω_μ_list    = df_list |> Map(df -> df.sol  |> Map(datum -> datum.Ω_μ)                   |> collect)  |> collect
        
r_f_result   = r_list |> Map(r  -> r[end])      |> collect
r_min_result = r_list |> Map(r  -> minimum(r))  |> collect
int_A_result = A_list |> Map(A -> ( sum(A) - (A[1] + A[end]) / 2 ) * 0.01 ) |> collect
t_f_result   = t_list |> Map(t -> t[end])       |> collect
avg_A_result = int_A_result ./ t_f_result

#  s_sim = Results(r_f_result, r_min_result, t_list, x_list, y_list, A_list, r_list, σ_M_list)
       
# OX
# f_OX = plot()
# for i in 1:length(s_BPNG_list)
#     if r_min_result[i] < 4 && r_f_result[i] < 4
#         plot!(f_OX, i * ones(2), ones(2), markershape=:circle, markercolor=:blue, label=:false)
#     else
#         plot!(f_OX, i * ones(2), ones(2), markershape=:xcross, markercolor=:red, label=:false)
#     end
# end
# display(f_OX)
        
# xy
f_2D = plot(legend_column = -1)
for i in 1:length(s_BPNG_list)
    fig_print(x_list[i], y_list[i], [], :false, "\$x ~[\\textrm{m}]\$", "\$y ~[\\textrm{m}]\$", f_2D; save_file=0)
end
fig_print([], [], "Traj2D_Case$(case)", [], [], [], f_2D; ar_val=:equal)
# display(f_2D)

# xyz
f_3D = plot()
for i in 1:length(s_BPNG_list)
    if case == 1
        label_input = label_string[i]
    else
        label_input = :false
    end
    equal_AR_3D(x_list[i], y_list[i], z_list[i], label_input, f_3D)
    # plot!(f_3D, x_list[i], y_list[i], z_list[i], label_string[i])
end
plot!(f_3D, xlabel="\$x ~[\\textrm{m}]\$", ylabel="\$y ~[\\textrm{m}]\$", zlabel="\$z ~[\\textrm{m}]\$", camera=camera_angle)
fig_print([], [], "Traj_Case$(case)", [], [], [], f_3D; lgnd_val = :outertopright, lfs_val = 12, ar_val=:equal)

# display(f_3D)

# r
f_r = plot()
for i in 1:length(s_BPNG_list)
    if case in [1,2,3,6]
        label_input = label_string[i]
    else
        label_input = :false
    end
    fig_print(t_list[i], r_list[i], [], label_input, "\$t ~[\\textrm{s}]\$", "\$r ~[\\textrm{m}]\$", f_r; save_file=0)
end
fig_print([], [], "r_Case$(case)", [], [], [], f_r; lgnd_val = :bottomleft, lfs_val = 12)
# display(f_r)

# σ_M
f_σ = plot()
for i in 1:length(s_BPNG_list)
    fig_print(t_list[i], rad2deg.(σ_M_list[i]), [], :false, "\$t ~[\\textrm{s}]\$", "\$\\sigma ~[\\textrm{deg}]\$", f_σ; save_file=0)
    
    if case == 3
        plot!(f_σ, t_list[i][[1,end],1], rad2deg.(s_BPNG_list[i].σ_M_lim * ones(2)), label=:false,linewidth=0.5, linestyle=:dash, linecolor=:red)
    end
end
fig_print([], [], "sigma_Case$(case)", [], [], [], f_σ)
# display(f_σ)

# e
f_e = plot()
for i in 1:length(s_BPNG_list)
    if case in [1,6]
        label_input = label_string[i]
    else
        label_input = :false
    end
    if case == 6 && isequal(s_BPNG_list[i].Bias, Bias_Quaternion_IACG)
        e_label = "\$\\varepsilon ~[\\textrm{deg}]\$"
    else
        e_label = "\$e ~[\\textrm{deg}]\$"
    end
    fig_print(t_list[i], rad2deg.(e_v̂_f_list[i]), [], label_input, "\$t ~[\\textrm{s}]\$", e_label, f_e; save_file=0)
end
fig_print([], [], "e_Case$(case)", [], [], [], f_e; lgnd_val = :topright)
# display(f_e)

# A = ||a_M||
f_A = plot()
for i in 1:length(s_BPNG_list)
    fig_print(t_list[i], A_list[i], [], :false, "\$t ~[\\textrm{s}]\$", "\$\\left||\\mathbf{a}\\right|| ~[\\textrm{m/s}^{2}]\$", f_A; save_file=0)
end
plot!(f_A, legend=:topleft)
fig_print([], [], "A_Case$(case)", [], [], [], f_A)
# display(f_A)

# A_bias = ||a_M_bias||
f_A_bias = plot()
for i in 1:length(s_BPNG_list)
    fig_print(t_list[i], A_bias_list[i], [], :false, "\$t ~[\\textrm{s}]\$", "\$\\left||\\mathbf{a}_{bias}\\right|| ~[\\textrm{m/s}^{2}]\$", f_A_bias; save_file=0)
end
fig_print([], [], "A_bias_Case$(case)", [], [], [], f_A_bias)
# display(f_A_bias)

# Ω_bias = ||ω_bias||
f_Ω_bias = plot()
for i in 1:length(s_BPNG_list)
    if case in 4:5
        label_input = label_string[i]
    else
        label_input = :false
    end
    fig_print(t_list[i], Ω_bias_list[i], [], label_input, "\$t ~[\\textrm{s}]\$", "\$\\left||\\mathbf{\\omega}_{bias}\\right|| ~[\\textrm{rad/s}]\$", f_Ω_bias; save_file=0)
end
fig_print([], [], "omega_bias_Case$(case)", [], [], [], f_Ω_bias; lgnd_val = :topright, lfs_val = 12)
# display(f_Ω_bias)

# μ
if case == 5
    f_μ = plot()
    for i in 1:length(s_BPNG_list)
        fig_print(t_list[i], rad2deg.(μ_list[i]), [], :false, "\$t ~[\\textrm{s}]\$", "\$\\mu ~[\\textrm{deg}]\$", f_μ; save_file=0)
    end
    fig_print([], [], "mu_Case$(case)", [], [], [], f_μ)
    # display(f_μ)
end

# Ω_μ
if case in [4,5]
    f_Ω_μ = plot()
    for i in 1:length(s_BPNG_list)
        fig_print(t_list[i], Ω_μ_list[i], [], :false, "\$t ~[\\textrm{s}]\$", "\$\\Omega_{\\mu} ~[\\textrm{rad/s}]\$", f_Ω_μ; save_file=0)
    end
    fig_print([], [], "Omega_mu_Case$(case)", [], [], [], f_Ω_μ)
    # display(f_Ω_μ)
end

@show int_A_result;
@show avg_A_result;

println("EOS")
