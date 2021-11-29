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
begin
	include("sim_missiles.jl")

	# To avoid redefinition of a named (or an anonymous) function at each time tuning parameters and running simulation,
	# place the custom definition of a bias command function in a different cell from the cell containing sim_plot.

	custom_bias(s,x,t) = zeros(3)                           # option 1
	# function custom_bias(s,x,t)  return zeros(3) end      # option 2
	# custom_bias = (s,x,t) -> zeros(3)                     # option 3
end

# ╔═╡ b4aa1cc3-07de-4717-a39e-8a9c4b1c8512
md"""
# Missile Guidance Simulation
2D/3D Vector-Form Biased Proportional Navigation Guidance Law is implemented.
The bias command currently implemented considers target being stationary.

## Vector BPNG Equations
* Three-dimensional lead angle
$\sigma = \cos^{-1}\left(\hat{\mathbf{r}} \cdot \hat{\mathbf{v}}\right) \in \left[0,\pi\right]$

* Manoeuvre plane unit normal vector
$\hat{\mathbf{k}} = \frac{\mathbf{r} \times \mathbf{v}}{\left\|\mathbf{r} \times \mathbf{v}\right\|}$

* Predicted final velocity direction vector
$\begin{aligned}
	\hat{\mathbf{v}}_{f_{pred}} &= \hat{\mathbf{v}} \cos\left(\frac{N}{N-1}\sigma\right) - \left(\hat{\mathbf{k}}\times\hat{\mathbf{v}}\right)\sin\left(\frac{N}{N-1}\sigma\right)\\
	&= \hat{\mathbf{r}} \cos\left(\frac{1}{N-1}\sigma\right) - \left(\hat{\mathbf{k}}\times\hat{\mathbf{r}}\right)\sin\left(\frac{1}{N-1}\sigma\right)
\end{aligned}$

* Predicted final velocity direction error angle
$e = \cos^{-1}\left(\hat{\mathbf{v}}_{f_{pred}} \cdot \hat{\mathbf{v}}_{f_{d}}\right) \in \left[0,\pi\right]$

* Desired angular velocity of predicted final velocity direction vector
$\boldsymbol{\omega}_{f_{d}} = -K_{e}\left(r\right) \dot{r} e   \frac{\hat{\mathbf{v}}_{f_{pred}} \times \hat{\mathbf{v}}_{f_{d}}}{\left\|\hat{\mathbf{v}}_{f_{pred}} \times \hat{\mathbf{v}}_{f_{d}}\right\|} + K_{roll}\left(r\right) \hat{\mathbf{v}}_{f_{d}}$

* Bias angular velocity command
$\boldsymbol{\omega}_{bias} = K_{\sigma}\left(\sigma\right) \left[-\left(N-1\right) \left(\boldsymbol{\omega}_{f_{d}}\cdot \hat{\mathbf{k}}\right) \hat{\mathbf{k}}  + \sin\sigma \,\boldsymbol{\omega}_{f_{d}} \cdot \left\{\hat{\mathbf{r}} + \cot\left(\frac{1}{N-1}\sigma\right)\left(\hat{\mathbf{k}}\times\hat{\mathbf{r}}\right)\right\} \left(\hat{\mathbf{v}} \times \hat{\mathbf{k}}\right)\right]$

* LOS angular velocity
$\boldsymbol{\omega}_{r}^{r\perp} = \frac{\mathbf{v} \times \mathbf{r}}{\left\| \mathbf{r} \right\|^{2}} = -\frac{v\sin\sigma}{r}\hat{\mathbf{k}}$

* Lateral acceleration command
$\mathbf{a}_{BPNG}^{v\perp} = \left(N\boldsymbol{\omega}_{r}^{r\perp} + \boldsymbol{\omega}_{bias}\right) \times \mathbf{v}$

"""

# ╔═╡ 9fcf0c9b-d559-402f-883c-c53f61815c94
begin
	md"""
	 $\alpha$ = $(@bind α Slider(0.01:0.01:1; default=1, show_value=true))
	
	 $n$ = $(@bind n Slider(0 : 0.1 : 10; default=1, show_value=true))
	
	 $\gamma_{f_{d}}$ [deg] = $(@bind γ_f_d_deg Slider(-90:1:90; default=0, show_value=true))
	
	 $\chi_{f_{d}}$ [deg] = $(@bind χ_f_d_deg Slider(0:1:180; default=180, show_value=true))
	
	 $\sigma_{M_{\lim}}$ [deg] = $(@bind σ_M_lim_deg Slider(0:1:90; default=60, show_value=true))
	"""
end

# ╔═╡ b75636f4-b7f2-4356-9e2a-d90942c0071e
begin
	s_BPNG.Bias    = Bias_IACG_StationaryTarget
	s_BPNG.s_Bias  = ComponentArray(α = 1, δ = 0.01, n = n, r_ref = 10E3, k = 9, m = 10)
	s_BPNG.σ_M_lim = deg2rad(σ_M_lim_deg)
	γ_f_d = deg2rad(γ_f_d_deg)
	χ_f_d = deg2rad(χ_f_d_deg)
	s_BPNG.v̂_f_d = [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)]
	df, f = sim_plot(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG, A_M_max)
	f
end

# ╔═╡ bce435b0-c6d4-4ac7-93aa-727ed08c73e5
begin
	t = df.time
	r = df.sol |> Map(datum -> datum.r) |> collect
	plot(t,r, label = :false)
	annotate!((0.9,0.9), "min r = "*string(minimum(r)) )
end

# ╔═╡ Cell order:
# ╟─b4aa1cc3-07de-4717-a39e-8a9c4b1c8512
# ╟─02329bbc-15a4-11ec-2522-df175bfb53e0
# ╠═50becda1-a0e6-4e4c-98ca-a72992a5fc42
# ╠═b9c304e7-9034-44d4-a4a6-a5d34d0baf58
# ╟─9fcf0c9b-d559-402f-883c-c53f61815c94
# ╠═b75636f4-b7f2-4356-9e2a-d90942c0071e
# ╟─bce435b0-c6d4-4ac7-93aa-727ed08c73e5
