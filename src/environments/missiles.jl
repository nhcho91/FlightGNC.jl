abstract type AbstractMissile <: AbstractEnv end

#-----------------------------------------
# Single Point Mass Vehicle Kinematics
#-----------------------------------------
struct PointMassMissile <: AbstractMissile
    dim::Int  # @assert dim == 2 || dim == 3
    A_max::Float64
end

function State(env::PointMassMissile)
    @unpack dim = env
    return function (p = zeros(dim), v = zeros(dim))
        ComponentArray(p = p, v = v)
    end
end

function Dynamics!(env::PointMassMissile)
    @unpack A_max = env
    @Loggable function dynamics!(dx, x, params, t; u)
        @unpack p, v = x
        @log p
        @log v
        dx.p = v
        if norm(u) >= A_max
            u = normalize(u) * clamp(norm(u), 0, A_max)
        end
        dx.v = u
        @log u
    end
end


#-----------------------------------------
# 1 Pursuer on 1 Evador Engagement
#-----------------------------------------
struct PursuerEvadorMissile <: AbstractMissile
    pursuer::PointMassMissile
    evador::PointMassMissile
end

function State(env::PursuerEvadorMissile)
    @unpack pursuer, evador = env
    return function (x0_pursuer, x0_evador)
        ComponentArray(pursuer = x0_pursuer, evador = x0_evador)
    end
end

function Dynamics!(env::PursuerEvadorMissile, s_guidance::BPNG)
    @unpack pursuer, evador = env
    @Loggable function dynamics!(dx, x, params, t; u_pursuer, u_evador)
        @onlylog r = norm(x.pursuer.p - x.evador.p)
        @onlylog σ_M = acos(clamp(dot(normalize(x.evador.p - x.pursuer.p), normalize(x.pursuer.v)), -1, 1))
        @onlylog a_M_bias, e_v̂_f, ω_bias, μ, Ω_μ = Bias_IACG_StationaryTarget(s_guidance, x, t)

        @nested_log :pursuer Dynamics!(pursuer)(dx.pursuer, x.pursuer, nothing, t; u = u_pursuer)
        @nested_log :evador Dynamics!(evador)(dx.evador, x.evador, nothing, t; u = u_evador)
    end
end


#-----------------------------------------
# Single Point Mass Vehicle Vertical Plane Dynamics
#-----------------------------------------
struct VerticalPlaneDynamicMissile <: AbstractMissile
    m::Float64
    S::Float64
    C_D_0::Float64
    K::Float64
    ρ₀::Float64
    H::Float64
    g::Float64
end

function State(env::VerticalPlaneDynamicMissile)
    return function (x₀, h₀, γ₀, V₀)
        ComponentArray(x = x₀, h = h₀, γ = γ₀, V = V₀)
    end
end

function Params(env::VerticalPlaneDynamicMissile)
    return function (R₀, c₀)
        ComponentArray(R = R₀, c = c₀)
    end
end

function Dynamics!(env::VerticalPlaneDynamicMissile)
    @unpack m, S, C_D_0, K, ρ₀, H, g = env
    @Loggable function dynamics!(ds, s, params, t; u)
        @unpack x, h, γ, V = s
        @unpack R, c = params
        @log x, h, γ, V, u
        @onlylog R, c

        ρ = ρ₀ * exp(-h / H)
        C_L = m * u / (1 / 2 * ρ * V^2 * S)
        C_D = C_D_0 + K * C_L^2

        ds.x = V * cos(γ)
        ds.h = V * sin(γ)
        ds.γ = u / V - g / V * cos(γ)
        ds.V = -C_D / C_L * u - g * sin(γ)
    end
end



#-----------------------------------------
# Single Rigid Body Vehicle Vertical Plane Dynamics
#-----------------------------------------
# [Disclaimer] Source of model parameters  
# @inproceedings{Mracek_1997,
# 	author = {Mracek, Curtis P. and Cloutier, James R.},
# 	title = {{Full Envelope Missile Longitudinal Autopilot Design Using the State-Dependent Riccati Equation Method}},
# 	booktitle = {AIAA Guidance, Navigation, and Control Conference},
# 	year = {1997},
# 	month = {August},
# 	address = {New Orleans, LA, USA},
# 	doi = {10.2514/6.1997-3767}
# }


Base.@kwdef struct VerticalPlane_RigidBody_Missile <: AbstractMissile
    # s_aerocoeff::Vector{Float64} = [-0.3, 19.373, -31.023, -9.717, -1.948, 40.44, -64.015, 2.922, -11.803]
    # s_airframe::Vector{Float64} = [204.02, 247.439, 0.0409, 0.2286, 150.0, 0.7]
    # s_environment::Vector{Float64} = [9.8, 1.225, 8435.0, 1.4, 286.0, 288.15, 0.0065]
    
    s_aerocoeff::Vector{Float64} = ComponentArray(aₐ = -0.3, aₙ = 19.373, bₙ = -31.023, cₙ = -9.717, dₙ = -1.948, aₘ = 40.44, bₘ = -64.015, cₘ = 2.922, dₘ = -11.803)
    s_airframe::Vector{Float64} = ComponentArray(m = 204.02, I_yy = 247.439, S = 0.0409, d = 0.2286, ωₐ = 150.0, ζₐ = 0.7)
    s_environment::Vector{Float64} = ComponentArray(g = 9.8, ρ₀ = 1.225, H = 8435.0, γₐ = 1.4, Rₐ = 286.0, T₀ = 288.15, λ = 0.0065)
end

function State(env::VerticalPlane_RigidBody_Missile)
    return function (h₀, V₀, α₀, q₀, θ₀, δ₀, δ̇₀)
        ComponentArray(h = h₀, V = V₀, α = α₀, q = q₀, θ = θ₀, δ = δ₀, δ̇ = δ̇₀)
    end
end

function Dynamics!(env::VerticalPlane_RigidBody_Missile)
    @unpack s_aerocoeff, s_airframe, s_environment = env
    (aₐ, aₙ, bₙ, cₙ, dₙ, aₘ, bₘ, cₘ, dₘ)     = s_aerocoeff
    (m, I_yy, S, d, ωₐ, ζₐ)                 = s_airframe
    (g, ρ₀, H, γₐ, Rₐ, T₀, λ)               = s_environment

    @Loggable function dynamics!(ds, s, params, t; u)
        (h, V, α, q, θ, δ, δ̇) = s
        δ_c = u
        @log h, V, α, q, θ, δ, δ̇, δ_c
        # @onlylog 

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

        ds.h = V * sin(γ)
        ds.V = Q * S / m * (C_N * sin(α) + C_A * cos(α)) - g * sin(γ)
        ds.α = α̇
        ds.q = Q * S * d / I_yy * C_M
        ds.θ = q
        ds.δ = δ̇
        ds.δ̇ = -ωₐ^2 * (δ - δ_c) - 2.0f0 * ζₐ * ωₐ * δ̇
    end
end