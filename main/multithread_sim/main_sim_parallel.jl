using FlightGNC
using UnPack, ComponentArrays, Transducers
using LinearAlgebra, DifferentialEquations  # for callbacks
using Plots

function main(p_M_0::Vector, v_M_0::Vector, p_T_0::Vector, v_T_0::Vector, s_guidance::BPNG; A_M_max = Inf)
    @unpack dim, σ_M_lim = s_guidance
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
    function condition_stop(σ_M_lim)
        return function (u, t, integrator)
            (p_M, v_M, p_T, v_T) = (u.pursuer.p, u.pursuer.v, u.evador.p , u.evador.v)
            r = norm(p_T-p_M)
            ṙ = dot(p_T-p_M, v_T-v_M) / r
            σ_M = acos( clamp( dot(normalize(p_T-p_M), normalize(v_M)), -1, 1 ) )

            (r < 0.5) || (r < 10 && ṙ >= 0) || (σ_M > σ_M_lim)
        end
    end
    affect!(integrator) = terminate!(integrator)  # See DiffEq.jl documentation
    cb_stop    = DiscreteCallback(condition_stop(σ_M_lim), affect!)
    cb = CallbackSet(cb_stop)  # useful for multiple callbacks
        
    # Execute Simulation
    # prob: DE problem, df: DataFrame
    simulator = Simulator(x0, 
                        apply_inputs(Dynamics!(env, s_guidance);
                                    u_pursuer = BPNG_cmd(s_guidance),
                                    u_evador = (x, params, t) -> zeros(evador.dim));
                        Problem = :ODE,
                        solver = Tsit5(),
                        tf = t_sim_f)

    # Non-interactive simulation: solve approach (automatically reinitialised)
    @time df = solve(simulator; 
                    savestep = Δt, 
                    callback = cb, 
                    reltol = 1e-8, 
                    abstol = 1e-8)

    # # Execute Simulation
    # # prob: DE problem, df: DataFrame		
    # @time prob, df = FSimBase.sim(
    #                      x0,  # initial condition
    #                      apply_inputs(Dynamics!(env);
    #                                   u_pursuer = BPNG_cmd(s_guidance),
    #                                   u_evador = (x, params, t) -> zeros(evador.dim));  # dynamics!; apply_inputs is exported from FS and is so useful for systems with inputs
    #                      tf=t_sim_f,
    #                      savestep=Δt,  # savestep is NOT simulation step
    #                      solver=Tsit5(),
    #                      callback=cb,
    #                      reltol=1e-8,
    #                      abstol=1e-8
    #                     )  # sim is exported from FS

	return df
end


# -----------------------------------------------------------------------------------------
# Simulation Input: Initial Condition / Design parameters
dim     = 2

p_M_0   = zeros(dim)

V_M_0   = 300
γ_M_0 	= deg2rad(0)
χ_M_0 	= deg2rad(60)

p_T_0   = [5E3; 0; 0]; #[10E3; 5E3; 5E3]
v_T_0   = zeros(dim) # 100*[-1; 0; 0]

γ_f_d   = 0
χ_f_d 	= deg2rad(180)

if dim == 2
    γ_M_0   = 0
    γ_f_d   = 0
end

v_M_0   = V_M_0*[cos(γ_M_0)*sin(χ_M_0); cos(γ_M_0)*cos(χ_M_0); sin(γ_M_0)]
v̂_f_d   = [cos(γ_f_d)*sin(χ_f_d); cos(γ_f_d)*cos(χ_f_d); sin(γ_f_d)]

if dim == 2
    (p_M_0, v_M_0, p_T_0, v_T_0) = (p_M_0, v_M_0, p_T_0, v_T_0) .|> x -> x[1:2]
end



σ_M_lim = deg2rad(60)
A_M_max = 300

N       = 3
s_Bias  = ComponentArray(α = 1, δ = 0.01, n = 0, r_ref = 5E3, k = 4, m = 10)

s_BPNG  = BPNG(N, dim, σ_M_lim, v̂_f_d, Bias_IACG_StationaryTarget_2D, s_Bias) 
# Bias options: Bias_zero, Bias_IACG_StationaryTarget, Bias_IACG_StationaryTarget_2D

function sim_plot(p_M_0::Vector, v_M_0::Vector, p_T_0::Vector, v_T_0::Vector, s_BPNG::BPNG, A_M_max::Number)
    # Execute simulation
    df = main(p_M_0, v_M_0, p_T_0, v_T_0, s_BPNG; A_M_max = A_M_max)

    ts = df.time
    p_Ms = df.sol |> Map(datum -> datum.pursuer.p) |> collect
    p_Ts = df.sol |> Map(datum -> datum.evador.p)  |> collect

    p_Ms = hcat(p_Ms...)'
    p_Ts = hcat(p_Ts...)'


    # Plotting
    legend_string = ["Missile" "Target"]
    
    if dim == 2
        f = fig_print([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], [], legend_string, "x [m]", "y [m]"; ar_val = :equal, save_file = 0)
        # plot!(f, xlims=(0, 6E3), ylims=(-3E3, 3E3))

        # f = fig_print(ts, [p_Ms[:,1] p_Ts[:,2]], "", ["x" "y"], "t [s]", "pos [m]"; save_file = 0)

    elseif dim == 3	
        f_2D= fig_print([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], [], legend_string, "x [m]", "y [m]"; ar_val = :equal, save_file = 0)
        # display(f_2D)

        # f_3D = plot([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], [p_Ms[:,3] p_Ts[:,3]], label = legend_string, 
        #         xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]", aspect_ratio = :equal,
        #         camera =(45,30)
        #         )
        # f_3D = equal_AR_3D([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], [p_Ms[:,3] p_Ts[:,3]], [], f_3D)
        
        f_3D = equal_AR_3D([p_Ms[:,1] p_Ts[:,1]], [p_Ms[:,2] p_Ts[:,2]], [p_Ms[:,3] p_Ts[:,3]], legend_string)
        plot!(f_3D, xlabel = "x [m]", ylabel = "y [m]", zlabel = "z [m]", camera = (45,30))

        f = plot(f_2D, f_3D, layout = (2,1), size = (600, 1200))
    end
    
    f = fig_print([], [], "Trajectory", [], [], [], f; ar_val = :equal) 	# To just save the result
    # display(f)

    return df, f
end