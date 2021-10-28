include("main_Approx_V_f_Pred.jl")

io = open("V_f.txt", "w");
plot_on = 1

##

for Δt_replan_update in [1.0, 2.0, 5.0, 10.0, 20.0]
    for Approx_mode in 0:2
        s_guidance.Approx_mode = Approx_mode
        s_guidance.n           = 1.5

## Series 1: Initial Planning Only Simulation
        Δt_plan_update = 300.0
        df_plansim, t_plansim, x_plansim, h_plansim, γ_plansim, V_plansim, A_L_plansim, R_plansim, c_plansim = 
main(x₀, h₀, γ₀, V₀, Δt_plan_update, s_guidance)


## Series 2: Online Replanning Update Simulation
        Δt_plan_update = Δt_replan_update
        df_replansim, t_replansim, x_replansim, h_replansim, γ_replansim, V_replansim, A_L_replansim, R_replansim, c_replansim = 
main(x₀, h₀, γ₀, V₀, Δt_plan_update, s_guidance)

## Series 3: Initial Planner Output
x_plan, h_plan, γ_plan, V_plan, A_L_plan, R_plan, c_plan = plan_path(x₀, h₀, γ₀, V₀, s_guidance)

## Plotting
        if plot_on == 1
            label_string_xplot = ["PlanSim" "RePlanSim" "Plan" ]
            label_string_tplot = ["PlanSim" "RePlanSim"]

            f_x_h   = fig_print([x_plansim, x_replansim, x_plan], [h_plansim, h_replansim, h_plan],         "x_h_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_xplot, "\$x ~[\\textrm{m}]\$", "\$h ~[\\textrm{m}]\$"; save_file = 0) # ar_val = :equal
            plot!(f_x_h, legend = :bottomleft)
            f_x_h = fig_print([], [], "x_h_$(round(Int, Δt_plan_update))_$(Approx_mode)", [], [], [], f_x_h) 
            f_x_γ   = fig_print([x_plansim, x_replansim, x_plan], [γ_plansim, γ_replansim, γ_plan],     "x_gamma_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_xplot, "\$x ~[\\textrm{m}]\$", "\$\\gamma ~[\\textrm{deg}]\$")
            f_x_V   = fig_print([x_plansim, x_replansim, x_plan], [V_plansim, V_replansim, V_plan],         "x_V_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_xplot, "\$x ~[\\textrm{m}]\$", "\$V ~[\\textrm{m}/\\textrm{s}]\$")
            f_x_A_L = fig_print([x_plansim, x_replansim, x_plan], [A_L_plansim, A_L_replansim, A_L_plan], "x_A_L_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_xplot, "\$x ~[\\textrm{m}]\$", "\$A_{L} ~[\\textrm{m}/\\textrm{s}^{2}]\$")

            f_t_x   = fig_print([t_plansim, t_replansim], [x_plansim, x_replansim],         "t_x_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_tplot, "\$t ~[\\textrm{s}]\$", "\$x ~[\\textrm{m}]\$")
            f_t_h   = fig_print([t_plansim, t_replansim], [h_plansim, h_replansim],         "t_h_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_tplot, "\$t ~[\\textrm{s}]\$", "\$h ~[\\textrm{m}]\$")
            f_t_γ   = fig_print([t_plansim, t_replansim], [γ_plansim, γ_replansim],     "t_gamma_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_tplot, "\$t ~[\\textrm{s}]\$", "\$\\gamma ~[\\textrm{deg}]\$")
            f_t_V   = fig_print([t_plansim, t_replansim], [V_plansim, V_replansim],         "t_V_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_tplot, "\$t ~[\\textrm{s}]\$", "\$V ~[\\textrm{m}/\\textrm{s}]\$")
            f_t_A_L = fig_print([t_plansim, t_replansim], [A_L_plansim, A_L_replansim],   "t_A_L_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_tplot, "\$t ~[\\textrm{s}]\$", "\$A_{L} ~[\\textrm{m}/\\textrm{s}^{2}]\$")
            f_t_R   = fig_print([t_plansim, t_replansim], [R_plansim, R_replansim],         "t_R_$(round(Int, Δt_plan_update))_$(Approx_mode)", label_string_tplot, "\$t ~[\\textrm{s}]\$", "\$R\$")

            f_t_c = plot()
            fig_print(t_plansim,   c_plansim,     [], "PlanSim",   "\$t ~[\\textrm{s}]\$", "\$\\mathbf{c}\$", f_t_c; save_file = 0)
            fig_print(t_replansim, c_replansim,   [], "RePlanSim", "\$t ~[\\textrm{s}]\$", "\$\\mathbf{c}\$", f_t_c; save_file = 0)
            f_t_c = fig_print([], [], "t_c_$(round(Int, Δt_plan_update))_$(Approx_mode)", [], [], [], f_t_c) 
        end

        write(io, "Δt_replan_update: $(round(Int, Δt_plan_update)) / Approx_mode: $(Approx_mode) \n")
        write(io, "x_f: plansim = $(x_plansim[end]), replansim = $(x_replansim[end]), plan = $(x_plan[end])\n")
        write(io, "h_f: plansim = $(h_plansim[end]), replansim = $(h_replansim[end]), plan = $(h_plan[end])\n")
        write(io, "γ_f: plansim = $(γ_plansim[end]), replansim = $(γ_replansim[end]), plan = $(γ_plan[end])\n")
        write(io, "V_f: plansim = $(V_plansim[end]), replansim = $(V_replansim[end]), plan = $(V_plan[end])\n\n\n")
    end
end
close(io)

## -----------------------------------------------------------------------------------------
# Test
# -----------------------------------------------------------------------------------------
# ##
# status, J_opt, c_opt = Path_Generator(s_guidance, param.R, x₀, h₀, γ₀)
# V̂_f = Speed_Predictor(s_guidance, x₀, V₀, c_opt)

# ##
# R_list = -5:1:5
# V̂_f_list = R_list |> Map(R -> begin
#     status, J_opt, c_opt = Path_Generator(s_guidance, R, x₀, h₀, γ₀)
#     V̂_f = Speed_Predictor(s_guidance, x₀, V₀, c_opt)
# end) |> tcollect

##
# s₀     = State(env_model)(x₀, h₀, γ₀, V₀)
# R₀     = 1.0
# c₀     = zeros(p)
# param₀ = Params(env_model)(R₀, c₀)

# s_guidance.Approx_mode = 0
# param_exact  = IASCG_Planner(s_guidance, s₀, param₀)
# s_guidance.Approx_mode = 1
# s_guidance.n = 1.5
# param_approx = IASCG_Planner(s_guidance, s₀, param₀)

# ##
# x = x₀:1:x_f

# order = p-1:-1:0
# h_exact = x |> Map(x -> h_ref * ((x / x_ref).^order)' * param_exact.c) |> collect
# h_exact=hcat(h_exact...)'
# h_approx = x |> Map(x -> h_ref * ((x / x_ref).^order)' * param_approx.c) |> collect
# h_approx=hcat(h_approx...)'

# f_x_h = fig_print(x, [h_exact h_approx], "x_h_plan_comparison", ["exact" "approx"], "\$x ~[\\textrm{m}]\$", "\$h ~[\\textrm{m}]\$"; save_file=1)
# # fig_print([], [], "x_h", [], [], [], f_x_h; ar_val=:equal, lfs_val=6)
