using FlightGNC
using UnPack, ComponentArrays, Transducers
using LinearAlgebra, DifferentialEquations
using Plots

function main(x₀::Float64, h₀::Float64, γ₀::Float64, V₀::Float64, Δt_plan_update::Float64, s_guidance::IASCG)
    @unpack p, env_model = s_guidance
    @unpack m, S, C_D_0, K, ρ₀, H, g = env_model

    # Environments
    env     = VerticalPlaneDynamicMissile(m, S, C_D_0, K, ρ₀, H, g)

    # Initial condition
    s₀      = State(env)(x₀, h₀, γ₀, V₀)
    R₀      = 1.0
    c₀      = zeros(p)
    params₀ = Params(env)(R₀, c₀)

    # Simulation parameters
    Δt_save        = 0.01
    t_sim_f        = 100.0
    
    # callbacks
    function guidance_update(s_guidance::IASCG)
        return function (integrator)
            # println("t: $(integrator.t)")
            # println(integrator.p)
            if abs(integrator.u[1] - s_guidance.x_f) > 1E3
                integrator.p = IASCG_Planner(s_guidance, integrator.u, integrator.p) 
            end
        end
    end
    cb_plan_update = PeriodicCallback(guidance_update(s_guidance), Δt_plan_update; initial_affect = true)
    
    function condition_stop(s_guidance::IASCG)
        @unpack x_f = s_guidance
        return function (u, t, integrator)
            u.x - x_f
        end
    end
    affect!(integrator) = terminate!(integrator)
    cb_stop = ContinuousCallback(condition_stop(s_guidance), affect!)
    
    cb = CallbackSet(cb_plan_update, cb_stop)

    # Execute Simulation
    # prob: DE problem, df: DataFrame
    simulator = Simulator(s₀,  # initial condition
                        apply_inputs(Dynamics!(env);
                                    u = IASCG_cmd(s_guidance)),
                        params₀;
                        Problem = :ODE,
                        solver = Tsit5(),
                        tf = t_sim_f)

    # Non-interactive simulation: solve approach (automatically reinitialised)
    @time df = solve(simulator; 
                    savestep = Δt_save, 
                    callback = cb, 
                    reltol = 1e-8, 
                    abstol = 1e-8)

    # # Execute Simulation
    # # prob: DE problem, df: DataFrame		
    # @time prob, df = FSimBase.sim(
    #         s₀,  # initial condition
    #         apply_inputs(Dynamics!(env);
    #                      u = IASCG_cmd(s_guidance)),
    #         params₀;  # dynamics!; apply_inputs is exported from FS and is so useful for systems with inputs
    #         tf = t_sim_f,
    #         savestep = Δt_save,  # savestep is NOT simulation step
    #         solver = Tsit5(),
    #         callback = cb,
    #         reltol = 1e-8,
    #         abstol = 1e-8
    #        ) 

    t_sim   = df.time
    x_sim   = df.sol |> Map(datum -> datum.x)  |> collect
    h_sim   = df.sol |> Map(datum -> datum.h)  |> collect
    γ_sim   = rad2deg.( df.sol |> Map(datum -> datum.γ)  |> collect )
    V_sim   = df.sol |> Map(datum -> datum.V)  |> collect
    A_L_sim = df.sol |> Map(datum -> datum.u)  |> collect
    R_sim   = df.sol |> Map(datum -> datum.R)  |> collect
    c_sim   = df.sol |> Map(datum -> datum.c)  |> collect

    c_sim  = hcat(c_sim...)'

    return df, t_sim, x_sim, h_sim, γ_sim, V_sim, A_L_sim, R_sim, c_sim
end

function plan_path(x₀::Float64, h₀::Float64, γ₀::Float64, V₀::Float64, s_guidance::IASCG)
    @unpack p, Δx_pred, env_model = s_guidance
    @unpack g = env_model

    s₀     = State(env_model)(x₀, h₀, γ₀, V₀)
    R₀     = 1.0
    c₀     = zeros(p)
    param₀ = Params(env_model)(R₀, c₀)
    param_plan = IASCG_Planner(s_guidance, s₀, param₀)

    x_plan = x₀:Δx_pred:x_f
    order = p-1:-1:0
    path_s = x_plan |> Map(x -> begin
                            ϕ  = (x / x_ref).^order
                            ϕ′ = order .* [ϕ[2:end]; 0] / x_ref
                            ϕ″ = order .* (order .- 1) .* [ϕ[3:end]; zeros(2)] / x_ref^2
                            # path_deriv = h_ref * [ϕ'; ϕ′'; ϕ″'] * param_plan[:c]
                            path = [h_ref * dot(ϕ, param_plan[:c]) atand(h_ref * dot(ϕ′, param_plan[:c])) h_ref * dot(ϕ″, param_plan[:c])]
                        end) |> tcollect
    path_s  = vcat(path_s...)
    h_plan  = path_s[:,1]
    γ_plan  = path_s[:,2]
    h″_plan = path_s[:,3]

    prob    = ODEProblem(Speed_Dynamics(s_guidance, param_plan[:c]), V₀, (x₀, x_f))
    sol     = DifferentialEquations.solve(prob, Tsit5(), reltol = 1e-6, abstol = 1e-6, saveat = x_plan)
    V_plan  = sol.u

    A_L_plan = ( (h″_plan .* ((cos.(deg2rad.(γ_plan))).^2) .* (V_plan.^2)) .+ g ) .* cos.(deg2rad.(γ_plan))

    return x_plan, h_plan, γ_plan, V_plan, A_L_plan, param_plan[:R], param_plan[:c]
end


# -----------------------------------------------------------------------------------------
# Simulation Input Default Values
# -----------------------------------------------------------------------------------------
# Algorithm parameters
p       = 10
x_ref   = 40E3
h_ref   = 3E3
h_des   = 1E3
Approx_mode = 1
n       = 1.5
Δx_pred = 1
x_f     = x_ref
h_f_d   = 0E3
γ_f_d   = deg2rad(-30)
V_f_d   = 400.0

# environmental parameters
m       = 544
S       = 0.258
C_D_0   = 0.126
K       = 0.370
ρ₀      = 1.225
H       = 8435
g       = 9.805

env_model  = VerticalPlaneDynamicMissile(m, S, C_D_0, K, ρ₀, H, g)

# initial conditions
R   = 1.0
x₀  = 0.0
h₀  = 3E3
γ₀  = deg2rad(5)
V₀  = 1000.0

s_guidance = IASCG(p, x_ref, h_ref, h_des, Approx_mode, n, x₀, V₀, Δx_pred, x_f, h_f_d, γ_f_d, V_f_d, env_model)
