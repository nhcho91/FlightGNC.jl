abstract type AbstractMissile <: AbstractEnv end

"""
Single Point Mass Vehicle Kinematics
"""
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
            u = normalize(u) * clamp(norm(u), -A_max, A_max)
        end
        dx.v = u
        @log u
    end
end


"""
1 Pursuer on 1 Evador Engagement
"""
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

function Dynamics!(env::PursuerEvadorMissile)
    @unpack pursuer, evador = env
    @Loggable function dynamics!(dx, x, params, t; u_pursuer, u_evador)
        @onlylog r   = norm(x.pursuer.p - x.evador.p)
        @onlylog σ_M = acos( clamp( dot( normalize(x.evador.p - x.pursuer.p), normalize(x.pursuer.v) ), -1, 1 ) )
        
        @nested_log :pursuer Dynamics!(pursuer)(dx.pursuer, x.pursuer, nothing, t; u = u_pursuer)
        @nested_log :evador Dynamics!(evador)(dx.evador, x.evador, nothing, t; u = u_evador)
    end
end
