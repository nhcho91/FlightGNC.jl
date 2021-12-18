using ContinuousTimePolicyGradients
using DiffEqFlux, ComponentArrays, LinearAlgebra, JLD2, OrdinaryDiffEq

function main()
    # model + problem parameters
    (aₐ, aₙ, bₙ, cₙ, dₙ, aₘ, bₘ, cₘ, dₘ) = Float32.([-0.3, 19.373, -31.023, -9.717, -1.948, 40.44, -64.015, 2.922, -11.803])
    (m, I_yy, S, d, ωₐ, ζₐ) = Float32.([204.02, 247.439, 0.0409, 0.2286, 150.0, 0.7])
    (g, ρ₀, H, γₐ, Rₐ, T₀, λ) = Float32.([9.8, 1.225, 8435.0, 1.4, 286.0, 288.15, 0.0065])
    (a_max, α_max, δ_max, q_max, M_max, h_max) = Float32.([100.0, deg2rad(30), deg2rad(25), deg2rad(5), 4, 11E3])

    (k_a, k_δ, k_R) = Float32.([10.0, 0.01, 0.0])

    # dynamic model
    dim_x = 7
    function dynamics_plant(t, x, u)
        (h, V, α, q, θ, δ, δ̇) = x
        δ_c = u[1]

        ρ = ρ₀ * exp(-h / H)
        Vₛ = sqrt(γₐ * Rₐ * (T₀ - λ * h))
        M = V / Vₛ
        γ = θ - α
        Q = 0.5f0 * ρ * V^2

        C_A = aₐ
        C_N = aₙ * α^3 + bₙ * α * abs(α) + cₙ * (2.0f0 - M / 3.0f0) * α + dₙ * δ
        C_M = aₘ * α^3 + bₘ * α * abs(α) + cₘ * (-7.0f0 + 8.0f0 * M / 3.0f0) * α + dₘ * δ

        α̇ = Q * S / m / V * (C_N * cos(α) - C_A * sin(α)) + g / V * cos(γ) + q
        a_z = V * (q - α̇)

        dx = [V * sin(γ);
            Q * S / m * (C_N * sin(α) + C_A * cos(α)) - g * sin(γ);
            α̇;
            Q * S * d / I_yy * C_M;
            q;
            δ̇;
            -ωₐ^2 * (δ - δ_c) - 2.0f0 * ζₐ * ωₐ * δ̇]
        y = [a_z]
        return dx, y
    end

    dim_x_c = 1
    function dynamics_controller(t, x_c, x, r, p_NN, policy_NN)
        (h, V, α, q, θ, δ, _) = x
        a_z_cmd = r[1]
        x_int   = x_c[1]

        ρ = ρ₀ * exp(-h / H)
        Vₛ = sqrt(γₐ * Rₐ * (T₀ - λ * h))
        M = V / Vₛ
        γ = θ - α
        Q = 0.5f0 * ρ * V^2

        C_A = aₐ
        C_N = aₙ * α^3 + bₙ * α * abs(α) + cₙ * (2.0f0 - M / 3.0f0) * α + dₙ * δ

        α̇ = Q * S / m / V * (C_N * cos(α) - C_A * sin(α)) + g / V * cos(γ) + q
        a_z = V * (q - α̇)

        K_A, K_I, K_R = policy_NN([α / α_max; M / M_max; h / h_max], p_NN)

        dx_c = [-K_A * (a_z - a_z_cmd) + q - a_z_cmd / V]
        u    = [-K_I * x_int - K_R * q]
        return dx_c, u
    end

    # cost definition
    function cost_running(t, x, y, u, r)
        a_z     = y[1]
        δ_c     = u[1]
        a_z_cmd = r[1]
        return k_a * ((a_z - a_z_cmd) / a_max)^2 + k_δ * (δ_c / δ_max)^2
    end

    function cost_terminal(x_f, r)
        a_z_cmd  = r[1]
        (_, y)   = dynamics_plant(0.0f0, x_f, Float32[0.0])
        a_z      = y[1]
        return ((a_z - a_z_cmd) / a_max)^2
    end

    function cost_regularisor(p_NN)
        return k_R * norm(p_NN)^2
    end

    # NN construction
    dim_NN_hidden = 64
    dim_NN_input  = 3
    dim_K = 3
    K_lb  = Float32.(-0.3*ones(3))
    K_ub  = zeros(Float32, 3)
    policy_NN = FastChain(
        FastDense(dim_NN_input, dim_NN_hidden, tanh),
        FastDense(dim_NN_hidden, dim_NN_hidden, tanh),
        FastDense(dim_NN_hidden, dim_K),
        (x, p) -> (K_ub - K_lb) .* σ.(x) .+ K_lb
    )

    # scenario definition
    ensemble_list = [ (Float32[h₀; V₀; zeros(5)], Float32[a_z_cmd])
                     for h₀      = 5E3:1E3:8E3
                     for V₀      = 7E2:1E2:9E2
                     for a_z_cmd = 0E1:2E1:1E2 ]
    t_span = Float32.((0.0, 3.0))

    scenario = ComponentArray(ensemble_list = ensemble_list, t_span = t_span, dim_x = dim_x, dim_x_c = dim_x_c)

    # NN training
    result, fwd_sol_nominal, loss_history = CTPG_train(dynamics_plant, dynamics_controller, cost_running, cost_terminal, cost_regularisor, policy_NN, scenario; sense_alg = InterpolatingAdjoint(autojacvec = ZygoteVJP()), ensemble_alg = EnsembleGPUArray(), saveat = 0.1f0)

    return result, fwd_sol_nominal, loss_history
end


@time result, fwd_sol_nominal, loss_history = main()
# jldsave("autopilot_saveat_0p1.jld2"; result, fwd_sol_nominal, loss_history)
