using FlightGNC
using UnPack, ComponentArrays, Transducers
using LinearAlgebra, DifferentialEquations  # for callbacks

function main(p_M_0::Vector, v_M_0::Vector, p_T_0::Vector, v_T_0::Vector, s_guidance::BPNG; A_M_max = Inf)
    @unpack dim = s_guidance
    # Environments
    pursuer = PointMassMissile(dim, A_M_max)
    evador = PointMassMissile(dim, Inf)
    env = PursuerEvadorMissile(pursuer, evador)
    
    # Initial condition
    x0_pursuer = State(pursuer)(p_M_0, v_M_0)
    x0_evador = State(evador)(p_T_0, v_T_0)
    x0 = State(env)(x0_pursuer, x0_evador)

    # Simulation parameters
    Δt      = 0.01
    t_sim_f = 50

    # callbacks
    function condition_stop(u, t, integrator)
        (p_M, v_M, p_T, v_T) = (u.pursuer.p, u.pursuer.v, u.evador.p , u.evador.v)
        r = norm(p_T-p_M)
        ṙ = dot(p_T-p_M, v_T-v_M) / r
        r < 0.1  || (r < 10 && ṙ >= 0)
    end
    affect!(integrator) = terminate!(integrator)  # See DiffEq.jl documentation
    cb_stop    = DiscreteCallback(condition_stop, affect!)
    cb = CallbackSet(cb_stop)  # useful for multiple callbacks
        
    # Execute Simulation
    # prob: DE problem, df: DataFrame		
    @time prob, df = FSimBase.sim(
                         x0,  # initial condition
                         apply_inputs(Dynamics!(env);
                                      u_pursuer = BPNG_cmd(s_guidance),
                                      u_evador = (x, params, t) -> zeros(evador.dim));  # dynamics!; apply_inputs is exported from FS and is so useful for systems with inputs
                         tf=t_sim_f,
                         savestep=Δt,  # savestep is NOT simulation step
                         solver=Tsit5(),
                         callback=cb,
                         reltol=1e-8,
                         abstol=1e-8
                        )  # sim is exported from FS

	return df
end

