using ModelingToolkit, DifferentialEquations, Plots
@variables t
D = Differential(t)

function dynamics(; name, τ = 0.01, a_max = 300.0, V₀ = 300.0, x_T = 5E3, y_T = 0.0)
    states = @variables x(t) y(t) γ(t) a(t) a_cmd(t) r(t) σ(t) λ(t) V(t)
    ps     = @parameters τ=τ a_max=a_max x_T=x_T y_T=y_T

    eqs = [
        D(x) ~ V*cos(γ)
        D(y) ~ V*sin(γ)
        D(γ) ~ a/V
        D(a) ~ (min(max(a_cmd,-a_max),a_max)-a)/τ
        r    ~ sqrt((x_T-x)^2 + (y_T-y)^2)
        λ    ~ atan(y_T-y, x_T-x)
        σ    ~ γ-λ
        V    ~ V₀
    ]

    return ODESystem(eqs, t, states, ps; name)
end

function guidance(; name, σ_lim = deg2rad(60), m = 10.0, δ = 0.01, k = 4.0, n = 1, p = 0.8, N=3.0, r₀ = 5E3, γ_f_d = -deg2rad(50))
    states = @variables r(t) γ(t) σ(t) λ(t) λ̇(t) V(t) e_γ_f(t) K_r(t) K_η(t) a_cmd(t)
    ps     = @parameters σ_lim=σ_lim m=m δ=δ k=k n=n p=p N=N r₀=r₀ γ_f_d=γ_f_d
    
    eqs = [
        λ̇       ~ -V*sin(σ)/r
        e_γ_f   ~ (N*λ-γ)/(N-1.0)-γ_f_d
        K_r     ~ (N-1.0+δ)*(1.0+r/r₀)^n
        K_η     ~ 1.0+k*(1.0-abs(sin(σ)/sin(σ_lim))^m)
        a_cmd   ~ N*V*λ̇ + V^2*cos(σ)*K_r*K_η/r*abs(e_γ_f)^p*sign(e_γ_f) 
    ]

    return ODESystem(eqs, t, states, ps; name) 
end

@named model    = dynamics()
@named command  = guidance() 

eqs_c = [
    model.V     ~ command.V
    model.r     ~ command.r
    model.γ     ~ command.γ
    model.σ     ~ command.σ
    model.λ     ~ command.λ
    model.a_cmd ~ command.a_cmd
]

sys = compose(ODESystem(eqs_c, name = :sys), model, command)
sys_simplified = structural_simplify(sys)

##
IC = [model.x => 0.0
      model.y => 0.0
      model.γ => deg2rad(55)
      model.a => 0.0]
t_span = (0.0, 30.0)
# pv = [command.n => 1
#       command.p => 1.0]
# prob = ODEProblem(sys_simplified, IC, t_span, pv; jac = true)

indexof(sym, syms) = findfirst(isequal(sym), syms)

function condition_stop(u,t,integrator)
    u_sys = states(sys_simplified)
    p_sys = parameters(sys_simplified)

    x     = u[indexof(model.x, u_sys)]
    y     = u[indexof(model.y, u_sys)]
    γ     = u[indexof(model.γ, u_sys)]
    x_T   = integrator.p[indexof(model.x_T, p_sys)]
    y_T   = integrator.p[indexof(model.y_T, p_sys)]
    σ_lim = integrator.p[indexof(command.σ_lim, p_sys)]

    r = sqrt((x_T-x)^2 + (y_T-y)^2)
    λ = atan(y_T-y, x_T-x)
    σ = γ-λ

    return  (r < 0.5) || (σ > σ_lim) #|| (r < 10 && ṙ >= 0)
end
affect!(integrator) = terminate!(integrator)
cb_stop = DiscreteCallback(condition_stop, affect!)
# @time sol = solve(prob, Tsit5(), callback = cb_stop)

## 
case = 2
if case == 1
    n_list = 1
    p_list = [0.6, 0.8, 1.0]
elseif case == 2
    n_list = [0, 1, 2]
    p_list = 0.9
end

pv_list = [[command.n => n, command.p => p] for n in n_list for p in p_list]

fs_val = 14
lw_val = 2.5
ls_list = [:solid :dash :dot :dashdot]

f_xy = plot()
f_r  = plot()
f_σ  = plot()
f_γ  = plot()
f_e_γ_f = plot()
f_a  = plot()
f_a_cmd  = plot()

for pv in pv_list
    prob = ODEProblem(sys_simplified, IC, t_span, pv; jac = true)
    @time sol = solve(prob, solver = Tsit5(), callback = cb_stop, reltol = 1e-8, abstol = 1e-8) 

    lgnd = "\$n = $(pv[1][2])\$, \$p = $(pv[2][2])\$"

    plot!(f_xy, sol, vars = (model.x, model.y), aspect_ratio = :equal, label = lgnd, xlabel = "\$x\$ [m]", ylabel = "\$y\$ [m]", linewidth = lw_val, linestyle = ls_list[findfirst(isequal(pv), pv_list)], legendfontsize = fs_val, guidefontsize = fs_val, xtickfontsize = fs_val, ytickfontsize = fs_val)

    plot!(f_r, sol, vars = model.r, label = lgnd, xlabel = "\$t\$ [s]", ylabel = "\$r\$ [m]", linewidth = lw_val, linestyle = ls_list[findfirst(isequal(pv), pv_list)], legendfontsize = fs_val, guidefontsize = fs_val, xtickfontsize = fs_val, ytickfontsize = fs_val)

    plot!(f_σ, sol, vars = rad2deg(model.σ), label = lgnd, xlabel = "\$t\$ [s]", ylabel = "\$\\sigma\$ [deg]", linewidth = lw_val, linestyle = ls_list[findfirst(isequal(pv), pv_list)], legendfontsize = fs_val, guidefontsize = fs_val, xtickfontsize = fs_val, ytickfontsize = fs_val, ylims = (0,60))

    plot!(f_γ, sol, vars = rad2deg(model.γ), label = lgnd, xlabel = "\$t\$ [s]", ylabel = "\$\\gamma\$ [deg]", linewidth = lw_val, linestyle = ls_list[findfirst(isequal(pv), pv_list)], legendfontsize = fs_val, guidefontsize = fs_val, xtickfontsize = fs_val, ytickfontsize = fs_val)

    plot!(f_e_γ_f, sol, vars = rad2deg(command.e_γ_f), label = lgnd, xlabel = "\$t\$ [s]", ylabel = "\$e_{\\gamma_{f}}\$ [deg]", linewidth = lw_val, linestyle = ls_list[findfirst(isequal(pv), pv_list)], legendfontsize = fs_val, guidefontsize = fs_val, xtickfontsize = fs_val, ytickfontsize = fs_val)

    plot!(f_a, sol, vars = model.a, label = lgnd, xlabel = "\$t\$ [s]", ylabel = "\$a\$ [m/s\$^{2}\$]", linewidth = lw_val, linestyle = ls_list[findfirst(isequal(pv), pv_list)], legendfontsize = fs_val, guidefontsize = fs_val, xtickfontsize = fs_val, ytickfontsize = fs_val)

    plot!(f_a_cmd, sol, vars = model.a_cmd, label = lgnd, xlabel = "\$t\$ [s]", ylabel = "\$a_{cmd}\$ [m/s\$^{2}\$]", linewidth = lw_val, linestyle = ls_list[findfirst(isequal(pv), pv_list)], legendfontsize = fs_val, guidefontsize = fs_val, xtickfontsize = fs_val, ytickfontsize = fs_val)
end

fig_dir = "Simulation/Julia_BPNG/Figures"
mkpath(fig_dir)

savefig(f_xy, joinpath(fig_dir, "Fig_xy_$(case).pdf"))
savefig(f_r, joinpath(fig_dir, "Fig_r_$(case).pdf"))
savefig(f_σ, joinpath(fig_dir, "Fig_sigma_$(case).pdf"))
savefig(f_γ, joinpath(fig_dir, "Fig_gamma_$(case).pdf"))
savefig(f_e_γ_f, joinpath(fig_dir, "Fig_e_gamma_f_$(case).pdf"))
savefig(f_a, joinpath(fig_dir, "Fig_a_$(case).pdf"))
savefig(f_a_cmd, joinpath(fig_dir, "Fig_a_cmd_$(case).pdf"))

display(f_xy)
display(f_r)
display(f_σ)
display(f_γ)
display(f_e_γ_f)
display(f_a)
display(f_a_cmd)