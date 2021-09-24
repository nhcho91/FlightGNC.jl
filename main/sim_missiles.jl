## Simulation Input: Initial Condition / Design parameters
using Plots
α       = 1
dim     = 2
N       = 3
s_Bias  = ComponentArray(δ = 0.01, n = 1, r_ref = 10E3, k = 9, m = 10)

V_M_0   = 300
γ_M_0 	= deg2rad(45)
χ_M_0_deg = 30
χ_f_d_deg = 180
σ_M_lim_deg = 60
χ_M_0 	= deg2rad(χ_M_0_deg)
χ_f_d 	= deg2rad(χ_f_d_deg)
σ_M_lim = deg2rad(σ_M_lim_deg)

function sim_plot(N, dim, α, χ_M_0, χ_f_d, σ_M_lim, s_Bias)
    if dim == 2
        (p_M_0, v_M_0)  = (zeros(dim), V_M_0*[sin(χ_M_0); cos(χ_M_0)])
        (p_T_0, v_T_0)  = ([5E3; 1E3], zeros(dim))
        γ_f_d = 0
    elseif dim == 3
        (p_M_0, v_M_0)  = (zeros(dim), V_M_0*[cos(γ_M_0)*sin(χ_M_0); cos(γ_M_0)*cos(χ_M_0); sin(γ_M_0)])
        (p_T_0, v_T_0)  = ([10E3; 5E3; 5E3], 100*[-1; 0; 0])
        γ_f_d = -0.01
    end
    v̂_f_d   = [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)]

    s_BPNG = BPNG(N, dim, α, σ_M_lim, v̂_f_d, Bias_zero, s_Bias) 
    # Bias = Bias_zero, Bias_IACG_StationaryTarget, Bias_IACG_StationaryTarget_2D

    # Execute simulation
    df = main(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max = 100)

    ts = df.time
    p_Ms = df.sol |> Map(datum -> datum.pursuer.p) |> collect
    p_Ts = df.sol |> Map(datum -> datum.evador.p)  |> collect

    p_Ms = hcat(p_Ms...)'
    p_Ts = hcat(p_Ts...)';


    ## Plotting
    legend_string = ["Missile" "Target"]
        
    f = fig_print([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], "", legend_string, "x [m]", "y [m]"; ar_val = :equal, save_file = 0, N_markers = 10)
    plot!(f, xlims=(0, 6E3), ylims=(-3E3, 3E3))
    f = fig_print([], [], "filename"; fig_handle = f, ar_val = :equal) 	# To just save the result
    display(f) 

    if dim == 3	
        f = plot([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], [p_Ms[:,3] p_Ts[:,3]], label = legend_string)
        f = fig_print([], [], "filename"; fig_handle = f, ar_val = :equal)
    end
end

# Parallel execution
# Threads.nthreads()

# n   = 100 # number of scenarios
# Ns  = 1:n |> Map(i -> rand(1)[1]) |> collect
# Ns= 1:0.1:3
# sim_results = Ns |> Map(N -> main(N)) |> tcollect