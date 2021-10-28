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
        @onlylog r   = norm(x.pursuer.p - x.evador.p)
        @onlylog σ_M = acos( clamp( dot( normalize(x.evador.p - x.pursuer.p), normalize(x.pursuer.v) ), -1, 1 ) )
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

        ρ   = ρ₀ * exp( -h / H )
        C_L = m * u / ( 1/2 * ρ * V^2 * S )
        C_D = C_D_0 + K * C_L^2

        ds.x = V * cos(γ)
        ds.h = V * sin(γ)
        ds.γ = u / V - g / V * cos(γ)
        ds.V = -C_D / C_L * u - g * sin(γ)
    end
end