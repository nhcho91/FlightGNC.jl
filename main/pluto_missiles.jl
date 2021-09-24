### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 02329bbc-15a4-11ec-2522-df175bfb53e0
begin
    import Pkg
    # careful: this is _not_ a reproducible environment
    # activate the global environment
    Pkg.activate()

    using PlutoUI
	using FlightGNC
	using ComponentArrays, Transducers
	using Plots
end

# ╔═╡ 50becda1-a0e6-4e4c-98ca-a72992a5fc42
include("main_missiles.jl")

# ╔═╡ b9c304e7-9034-44d4-a4a6-a5d34d0baf58
include("sim_missiles.jl")

# ╔═╡ b4aa1cc3-07de-4717-a39e-8a9c4b1c8512
md"""
# Missile Guidance Simulation
2D/3D Biased Proportional Navigation Guidance Law is implemented.

In 2D, lateral acceleration command is given as follows:

$\displaystyle a_{cmd} = NV\dot{\lambda} - \frac{K_{r}\left(r\right)K_{\eta}\left(\eta\right)}{r}e_{\gamma_{f}}$
"""

# ╔═╡ 9fcf0c9b-d559-402f-883c-c53f61815c94
begin
	md"""
	 $\alpha$ = $(@bind α Slider(0.01:0.01:1; default=1, show_value=true))
	
	 $n$ = $(@bind n Slider(0 : 0.1 : 10; default=1, show_value=true))
	
	 $\chi_{M_{0}}$ [deg] = $(@bind χ_M_0_deg Slider(0:1:180; default=30, show_value=true))
	
	 $\chi_{f_{d}}$ [deg] = $(@bind χ_f_d_deg Slider(0:1:180; default=180, show_value=true))
	
	 $\sigma_{M_{\lim}}$ [deg] = $(@bind σ_M_lim_deg Slider(0:1:90; default=60, show_value=true))
	"""
end

# ╔═╡ b75636f4-b7f2-4356-9e2a-d90942c0071e
begin
	s_Bias  = ComponentArray(δ = 0.01, n = n, r_ref = 10E3, k = 9, m = 10)
	sim_plot(3, 2, α, deg2rad(χ_M_0_deg), deg2rad(χ_f_d_deg), deg2rad(σ_M_lim), s_Bias)
end

# ╔═╡ bce435b0-c6d4-4ac7-93aa-727ed08c73e5
# begin
# 	# Initial Conditions / Design parameters
#     dim     = 3
#     N       = 3
#     s_Bias  = ComponentArray(δ = 0.01, n = n, r_ref = 10E3, k = 9, m = 10)
    
#     V_M_0 	= 300
# 	γ_M_0 	= deg2rad(45)
# 	χ_M_0 	= deg2rad(χ_M_0_deg)
# 	χ_f_d 	= deg2rad(χ_f_d_deg)
# 	σ_M_lim = deg2rad(σ_M_lim_deg)
#     if dim == 2
#         (p_M_0, v_M_0)  = (zeros(dim), V_M_0*[sin(χ_M_0); cos(χ_M_0)])
#         (p_T_0, v_T_0)  = ([5E3; 1E3], zeros(dim))
#         γ_f_d = 0
#     elseif dim == 3
#         (p_M_0, v_M_0)  = (zeros(dim), V_M_0*[cos(γ_M_0)*sin(χ_M_0); cos(γ_M_0)*cos(χ_M_0); sin(γ_M_0)])
#         (p_T_0, v_T_0)  = ([10E3; 5E3; 5E3], 100*[-1; 0; 0])
# 		γ_f_d = -0.01
#     end
#     v̂_f_d   = [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)]

# 	s_BPNG = BPNG(N, dim, α, σ_M_lim, v̂_f_d, Bias_zero, s_Bias) 
# 	# Bias = Bias_zero, Bias_IACG_StationaryTarget, Bias_IACG_StationaryTarget_2D
	
# 	# Execute simulation
# 	df = main(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max = 100)
	
# 	ts = df.time
#     p_Ms = df.sol |> Map(datum -> datum.pursuer.p) |> collect
#     p_Ts = df.sol |> Map(datum -> datum.evador.p)  |> collect

#     p_Ms = hcat(p_Ms...)'
#     p_Ts = hcat(p_Ts...)'
	
# 	# Plotting
# 	legend_string = ["Missile" "Target"]
	
# 	f = fig_print([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], "", legend_string, "x [m]", "y [m]"; ar_val = :equal, save_file = 0, N_markers = 10)
# 	# plot!(f, xlims=(0, 6E3), ylims=(-3E3, 3E3))
# 	f = fig_print([], [], "filename"; fig_handle = f) 	# To just save the result
	
# 	if dim == 3	
# 		f = plot([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], [p_Ms[:,3] p_Ts[:,3]], label = legend_string)
# 		f = fig_print([], [], "filename"; fig_handle = f, ar_val = :equal)
# 	end
# end

# ╔═╡ Cell order:
# ╟─b4aa1cc3-07de-4717-a39e-8a9c4b1c8512
# ╟─02329bbc-15a4-11ec-2522-df175bfb53e0
# ╠═50becda1-a0e6-4e4c-98ca-a72992a5fc42
# ╠═b9c304e7-9034-44d4-a4a6-a5d34d0baf58
# ╟─9fcf0c9b-d559-402f-883c-c53f61815c94
# ╠═b75636f4-b7f2-4356-9e2a-d90942c0071e
# ╟─bce435b0-c6d4-4ac7-93aa-727ed08c73e5
