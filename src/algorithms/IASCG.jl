abstract type AbstractIASCG <: AbstractController end

mutable struct IASCG <: AbstractIASCG
    p::Int
    x_ref::Float64
    h_ref::Float64
    h_des::Float64
    Approx_mode::Int
    n::Float64
    x₀::Float64
    V₀::Float64
    Δx_pred::Float64
    x_f::Float64
    h_f_d::Float64
    γ_f_d::Float64
    V_f_d::Float64
    env_model
end

function Path_Generator(s_guidance::IASCG, R, x_t, h_t, γ_t)
    @unpack p, x_ref, h_ref, h_des, x_f, h_f_d, γ_f_d = s_guidance
    
    order = p-1:-1:0
    ϕ_t  = (x_t / x_ref).^order
    ϕ_f  = (x_f / x_ref).^order
    dϕ_t = order .* [ϕ_t[2:end]; 0]
    dϕ_f = order .* [ϕ_f[2:end]; 0]

    A = [ϕ_t'; ϕ_f'; dϕ_t'; dϕ_f']
    b = [h_t; h_f_d; x_ref * tan(γ_t); x_ref * tan(γ_f_d)] / h_ref

    C = zeros(p, p)
    for i in 1:p
        for j in 1:p
            C[i,j] = 1 / (2 * p - (i + j - 1))
        end
    end
    S = (x_f * ϕ_f * ϕ_f'  - x_t * ϕ_t * ϕ_t') / x_ref .* C
    H = [zeros(p-2, 2) diagm(order[1:end-2] .* (order[1:end-2] .- 1));
         zeros(2,p)]
    P = S + (10.0)^R * H * S * H'
    P = Symmetric((P + P') / 2)
    q = -2 * h_des / h_ref * S[:,end]

    c = Convex.Variable(p)
    problem = minimize(quadform(c, P; assume_psd = true) + q' * c, A * c == b)
    Convex.solve!(problem, Mosek.Optimizer; silent_solver = true)
    # c_opt = evaluate(c)

    return problem.status, problem.optval, c.value
end


function Speed_Dynamics(s_guidance::IASCG, c)
    @unpack p, x_ref, h_ref, n, x₀, V₀, x_f, V_f_d, Approx_mode, env_model = s_guidance
    @unpack m, S, C_D_0, K, ρ₀, H, g = env_model
    return function (V, ps, x)
        order = p-1:-1:0
        ϕ  = (x / x_ref).^order
        ϕ′ = order .* [ϕ[2:end]; 0] / x_ref
        ϕ″ = order .* (order .- 1) .* [ϕ[3:end]; zeros(2)] / x_ref^2
        (h, h′, h″) = h_ref * [ϕ'; ϕ′'; ϕ″'] * c
        # γ  = atan(h′)
        # γ′ = h″ / (1 + h′^2)
    
        ρ  = ρ₀ * exp(-h / H)
        if Approx_mode == 0
            V̂_guess = V
        elseif Approx_mode == 1
            V̂_guess = V_f_d - (V_f_d - V₀) * ((x_f - x) / (x_f - x₀))^n
        elseif Approx_mode == 2
            V̂_guess = Inf
        end
        C_L = 2 * m / (ρ * S * sqrt(1 + h′^2)) * (h″ / (1 + h′^2) + g / V̂_guess^2)
        C_D = C_D_0 + K * C_L^2

        dV = - ρ * S * C_D * sqrt(1 + h′^2) / ( 2 * m ) * V - g / V * h′

        return dV
    end
end

function Speed_Predictor(s_guidance::IASCG, x, V, c)
    @unpack x_f = s_guidance
    # @unpack m, S, C_D_0, K, ρ₀, H, g = env_model

    # ----------------------------------------------------------------------------
    # Numerical Integration - Solving ODE
    # ---------------------------------------------------------------------------- 
    prob = ODEProblem(Speed_Dynamics(s_guidance, c), V, (x, x_f))
    sol = OrdinaryDiffEq.solve(prob, Tsit5(), reltol = 1e-6, abstol = 1e-6, dense = false, save_on = false, save_everystep = false, save_start = false, save_end = true)
    V̂_f = sol[1]

    # ----------------------------------------------------------------------------
    # Numerical Integration - Calculation of Definite Integrals
    # ----------------------------------------------------------------------------
    # x_s = x:Δx_pred:x_f
        
    # order = p-1:-1:0
    # path_deriv_s = x_s |> Map(x -> begin
    #                         ϕ  = (x / x_ref).^order
    #                         ϕ′ = order .* [ϕ[2:end]; 0] / x_ref
    #                         ϕ″ = order .* (order .- 1) .* [ϕ[3:end]; zeros(2)] / x_ref^2
    #                         path_deriv = h_ref * [ϕ'; ϕ′'; ϕ″'] * c
    #                     end) |> tcollect    
    # path_deriv_s = hcat(path_deriv_s...)'
    # (h_s, h′_s, h″_s) = path_deriv_s[[1,2,3],:]    

    # ρ_s = ρ₀ * exp.( - h_s / H)

    # a = -2 * (γ′_s .+ g ./ V̂_guess_s.^2) .* (C_D_0 .+ K .* C_L_s.^2) ./ C_L_s
    # b = -g * tan.(γ_s)   
    # Φ = exp(quadgk(a, x, x_f, rtol=1e-3))

    # V̂_f = sqrt( Φ[1] * V^2 + 2 * ∫Φb)

    return V̂_f
end



function zeroprob_IASCG_Planner(s_guidance::IASCG, x, h, γ, V)
    @unpack V_f_d = s_guidance
    return function (R)
        e_v_f = Speed_Predictor(s_guidance, x, V, Path_Generator(s_guidance, R, x, h, γ)[3]) - V_f_d

        return e_v_f
    end
end

function minbndprob_IASCG_Planner(s_guidance::IASCG, x, h, γ, V)
    @unpack V_f_d = s_guidance
    return function (R)
        e_v_f = Speed_Predictor(s_guidance, x, V, Path_Generator(s_guidance, R, x, h, γ)[3]) - V_f_d

        return e_v_f^2
    end
end

function IASCG_Planner(s_guidance::IASCG, s, params)
    @unpack p, V_f_d   = s_guidance
    @unpack x, h, γ, V = s
    @unpack R, c       = params
    
    # previous parameters -------------------------------------------
    # R_opt = R
    # c_opt = c

    # Root finding --------------------------------------------------
    # if c == zeros(p)
    #     R_opt = find_zero(zeroprob_IASCG_Planner(s_guidance, x, h, γ, V), (-8.0, 5.0), Bisection())
    # else
    #     R_opt = find_zero(zeroprob_IASCG_Planner(s_guidance, x, h, γ, V), R, Order1())
    # end

    # Minimisation over bounded interval ----------------------------
    res     = Optim.optimize(minbndprob_IASCG_Planner(s_guidance, x, h, γ, V), -8.0, 5.0)
    R_opt   = Optim.minimizer(res)

    c_opt = Path_Generator(s_guidance, R_opt, x, h, γ)[3]
    # V̂_f   = Speed_Predictor(s_guidance, x, V, c_opt)

    return ComponentArray(R = R_opt, c = vec(c_opt))
end

function IASCG_cmd(s_guidance::IASCG)
    @unpack p, x_ref, h_ref, env_model = s_guidance
    @unpack g = env_model
    return function (s, params, t)
        @unpack x, γ, V = s
        @unpack c = params

        order = p-1:-1:0
        ϕ  = (x / x_ref).^order
        # ϕ′ = order .* [ϕ[2:end]; 0] / x_ref
        ϕ″ = order .* (order .- 1) .* [ϕ[3:end]; zeros(2)] / x_ref^2
        h″ = h_ref * ϕ″' * c
        A_L_cmd = ( h″ * (cos(γ))^2 * V^2 + g ) * cos(γ)

        return A_L_cmd
    end
end

