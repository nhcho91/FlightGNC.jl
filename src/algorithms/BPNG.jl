abstract type AbstractBPNG <: AbstractController end

mutable struct BPNG <: AbstractBPNG
    N::Float64
    dim::Int
    σ_M_lim::Float64
    v̂_f_d::Vector
    Bias::Function
    s_Bias # parameters of bias command function
end

function BPNG_cmd(s_guidance::BPNG)
    @unpack N, dim, Bias = s_guidance
    return function (x, params, t)
        (p_M, v_M, p_T, v_T) = (x.pursuer.p, x.pursuer.v, x.evador.p, x.evador.v)
        if dim == 2
            (p_M, v_M, p_T, v_T) = vcat.((p_M, v_M, p_T, v_T), 0)
        end

        ω_r      = cross(p_T - p_M, v_T - v_M) / dot(p_T - p_M, p_T - p_M)
        a_M_PPNG = N * cross(ω_r, v_M)
        a_M_bias = Bias(s_guidance, x, t)[1]
        
        a_M      = a_M_PPNG + a_M_bias
        
        # if norm(a_M) >= A_M_max
        #     a_M = normalize(a_M) * clamp(norm(a_M), 0, A_M_max)
        # end
        if dim == 2
            a_M = a_M[1:2]
        end

        return a_M
    end
end

function Bias_IACG_StationaryTarget(s_guidance::BPNG, x, t)
    @unpack N, dim, σ_M_lim, v̂_f_d, s_Bias = s_guidance
    @unpack α, δ, n, r_ref, k, m, k̂_d, i_Ω_μ, i_σ_M_lim = s_Bias
    (p_M, v_M, p_T, v_T) = (x.pursuer.p, x.pursuer.v, x.evador.p, x.evador.v)
    if dim == 2
        (p_M, v_M, p_T, v_T) = vcat.((p_M, v_M, p_T, v_T), 0)
    end
    # @assert norm(v_T) < 0.1

    r        = norm(p_T - p_M)
    r̂        = normalize(p_T - p_M)
    ṙ        = dot(r̂, v_T - v_M)
    v̂_M      = normalize(v_M)
    σ_M      = atan(norm(cross(p_T - p_M, v_M)), dot(p_T - p_M, v_M))
    # σ_M      = acos( clamp( dot(r̂, v̂_M), -1, 1 ) )
    # atan-based implementation is more robustly accurate than using acos for computation
    k̂        = normalize(cross(p_T - p_M, v_M))

    v̂_f_pred = v̂_M * cos(N / (N - 1) * σ_M) - cross(k̂, v̂_M) * sin(N / (N - 1) * σ_M)
    # alternative: v̂_f_pred =   r̂ * cos(1 / (N - 1) * σ_M) -   cross(k̂, r̂) * sin(1 / (N - 1) * σ_M)
    e_v̂_f    = atan(norm(cross(v̂_f_pred, v̂_f_d)), dot(v̂_f_pred, v̂_f_d))

    # Design parameters for vec_BPNG
    k_e      = (n + 1) / r

    μ = 0
    if i_Ω_μ == 1
        # desired manoeuvre plane alignment
        k_μ  = (m + 1) / r
        μ    = atan(dot(cross(cross(v̂_f_pred, v̂_f_d), k̂_d), v̂_f_d), dot(cross(v̂_f_pred, v̂_f_d), k̂_d))
        if abs(μ) > pi / 2 
            μ    = atan(dot(cross(cross(v̂_f_pred, v̂_f_d), -k̂_d), v̂_f_d), dot(cross(v̂_f_pred, v̂_f_d), -k̂_d))
        end 
        Ω_μ  = k_μ * abs(ṙ) * μ

    elseif i_Ω_μ == 2
        # evasion: monotonic decrease
        Ω_μ  = δ * ( r / r_ref )^m
    elseif i_Ω_μ == 3
        # evasion: sinusoidal wave
        Ω_μ  = δ * sin(( r / r_ref )^m)
    elseif i_Ω_μ == 4
        # evasion: square wave
        Ω_μ  = δ * ifelse(mod(( r / r_ref )^m, 1) < θ, 1.0, -1.0)
    elseif i_Ω_μ == 5
        # evasion: sawtooth wave
        Ω_μ  = δ * rem(( r / r_ref )^m, 1.0, RoundNearest) * 2
    elseif i_Ω_μ == 6
        # evasion: triangular wave
        modr = mod(( r / r_ref )^m + 1/4, 1.0)
        Ω_μ  = δ * ifelse(modr < 1/2, 4*modr - 1, -4*modr + 3)
    elseif i_Ω_μ == 7
        # evasion: uniform random number (s_rand: a fixed sequence of m uniform randoms in (-1,1) generated externally)
        # Ω_μ  = δ * s_rand[1 + mod(floor(Int, t / r_ref), m)]
    else
        Ω_μ  = 0.0
    end

    if i_σ_M_lim == 1
        # lead angle limiter
        k_σ  = 2 / pi * acos( clamp((σ_M / σ_M_lim)^k , -1, 1) )
    else
        k_σ  = 1
    end

    ω_f      = k_e * abs(ṙ) * e_v̂_f^α * normalize(cross(v̂_f_pred, v̂_f_d)) + Ω_μ * v̂_f_d
    ω_bias   = k_σ * ( -(N - 1) * dot(ω_f, k̂) * k̂  + sin(σ_M) * dot(ω_f, r̂ + cot(1 / (N - 1) * σ_M) * cross(k̂, r̂)) * cross(v̂_M, k̂) )
    
    a_M_bias = cross(ω_bias, v_M)

    return a_M_bias, e_v̂_f, ω_bias, μ, Ω_μ
end



#--------------------------------------Legacy codes--------------------------------------
function Bias_IACG_StationaryTarget_2D(s_guidance::BPNG, x, t)
    @unpack N, dim, σ_M_lim, v̂_f_d, s_Bias = s_guidance
    @unpack α, δ, n, r_ref, k, m = s_Bias
    # @assert dim == 2
    (p_M, v_M, p_T, v_T) = (x.pursuer.p, x.pursuer.v, x.evador.p, x.evador.v)
    r       = norm(p_T - p_M)
    λ       = atan(p_T[2] - p_M[2], p_T[1] - p_M[1])   

    V_M     = norm(v_M)
    γ_M     = atan(v_M[2], v_M[1])
    σ_M     = γ_M - λ
    η       = sin(σ_M)
    η_lim   = sin(σ_M_lim)

    # V_T     = norm(v_T)
    # γ_T     = atan(v_T[2], v_T[1])
    # σ_T     = γ_T - λ

    γ_f_d   = atan(v̂_f_d[2], v̂_f_d[1])

    # λ̇       = (V_T*sin(σ_T) - V_M*sin(σ_M)) / r        
    e_γ_f   = γ_M - N / (N - 1) * σ_M - γ_f_d
    if α > 0.99
        e_γ_f_fbk = e_γ_f
    else
        e_γ_f_fbk = sign(e_γ_f) * abs(e_γ_f)^α
    end

    K_r     = (N - 1 + δ) * (1 + r / r_ref)^n 
    K_eta   = 1 + k * (1 - (abs(η / η_lim))^m)
    u_aug = - K_r * K_eta / r * e_γ_f_fbk

    A_M_bias = - u_aug * V_M^2 * cos(σ_M)
    a_M_bias = A_M_bias * [-sin(γ_M); cos(γ_M); 0]

    return a_M_bias

    # if dim == 2  # identical to the codes above, thus can be deleted
    #     Design parameters for LCG_e_0
    #     K_e@vec_BPNG = K_r@LCG_e_0 / (N - 1) / r
    #     kwargs: δ = 0.01, n = 0, r_ref = Inf, k = 0, m = 0
    #     K_e = (N - 1 + δ) / (N - 1) / r * (1 + r / r_ref)^n
    #     Ω_μ = 0    
    #     K_σ = 1 + k * (1 - (abs(sin(σ_M) / sin(σ_M_lim)))^m)
    #
    #     ê_z = [0; 0; 1]
    #     if cross(p_T - p_M, v_M)[3] >= 0
    #         k̂ = ê_z
    #     else
    #         k̂ = -ê_z
    #     end
    #     v̂_f_pred = v̂_M * cos(N / (N - 1) * σ_M) - cross(k̂, v̂_M) * sin(N / (N - 1) * σ_M)
    #     e_v̂_f_signed = atan(dot(v̂_f_pred, cross(ê_z, v̂_f_d)), dot(v̂_f_pred, v̂_f_d))
    #     # alternative: e_v̂_f_signed = atan(v̂_f_pred[2], v̂_f_pred[1]) - atan(v̂_f_d[2], v̂_f_d[1])
    #     ω_f     = K_e(r, N, δ, n, r_ref) * ṙ  * sign(e_v̂_f_signed) * abs(e_v̂_f_signed)^α * ê_z
    #     ω_bias  = K_σ(σ_M, σ_M_lim, k, m) * ( -(N - 1) * ω_f )
    # end
end

function Bias_zero(s_guidance::BPNG, x, t)
    return zeros(3)
end