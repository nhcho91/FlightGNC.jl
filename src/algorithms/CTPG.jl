"""
CTPG_train(dynamics_plant::Function, dynamics_controller::Function, cost_running::Function, cost_terminal::Function, cost_regularisor::Function, policy_NN, scenario; 
solve_alg = Tsit5(), sense_alg = InterpolatingAdjoint(autojacvec = ZygoteVJP()), ensemble_alg = EnsembleThreads(), opt_1 = ADAM(0.01), opt_2 = LBFGS(), maxiters_1 = 100, maxiters_2 = 100, progress_plot = true, solve_kwargs...)


`CTPG_train()` provides a high-level interface for optimisation of the neural networks inside an ODE-represented dynamics based on Continuous-Time Policy Gradient (CTPG) methods that belong to the adjoint sensitivity analysis techniques. The code implemented and the default values for keyword arguments are specified considering training of a neural controller as the main application. In the context herein, a neural controller refers to a dynamic controller that incorporates neural-network-represented components at some points in its mathematical description.

The code utilises the functionalities provided by the [DiffEqFlux.jl](https://github.com/SciML/DiffEqFlux.jl) and [DiffEqSensitivity.jl](https://github.com/SciML/DiffEqSensitivity.jl) packages, and the Automatic Differentiation (AD) capabilities provided by the [Zygote.jl](https://github.com/FluxML/Zygote.jl) package that is integrated in DiffEqFlux.jl. `CTPG_train()` presumes the consistency of the functions provided as its input arguments with the AD tool, hence, the dynamics and cost functions should maintain their transparence against AD tools.

The optimisation (training) problem minimises the cost function defined over deterministic samples of the initial plant state `x₀` and the reference `r` by performing ensemble simulation based on parallelised computation.

The signals are defined as described below:

- `t`: time
- `x`: plant state
- `y`: plant output (= sensor output)
- `x_c`: controller state
- `u`: plant input (= controller output)
- `r`: exogenous reference
- `x_aug`: augmented forward dynamics state (= `[x; x_c; ∫cost_running]`)
- `p_NN`: neural network parameter

The arguments should be provided as explained below:

- `dynamics_plant`: Describes the dynamics of the plant to be controlled. Input arguments `x` and `u` should be of Vector type.
- `dynamics_controller`: Describes the dynamics of the controller that includes neural networks components. Input arguments `x_c`, `y`, `r`, and `p_NN` should be of Vector type.
- `dynamics_sensor`: Describes the dynamics of the sensor that measures output variables fed to the controller. Input arguments `x` should be of Vector type: 
- `cost_running`: Describes the running cost defined as the integrand of the Lagrange-form continuous functional. Input arguments `x`, `y`, `u`, and `r` should be of Vector type.
- `cost_terminal`: Describes the terminal cost defined as the Mayer-form problem cost function. Defines a Bolza-form problem along with `cost_running`. Input arguments `x_f` and `r` should be of Vector type.
- `cost_regularisor`: Describes the regularisation term appended to the cost (loss) function. Input argument `p_NN` should be of Vector type.
- `policy_NN`: The neural networks entering into the controller dynamics. DiffEqFlux-based FastChain is recommended for its construction.
- `scenario`: Contains the parameters related with the ensemble-based training scenarios.
    - `ensemble`: A vector of the initial plant state `x₀` and the reference `r` constituting the trajectory realisations.
    - `t_span`: Time span for forward-pass integration
    - `t_save`: Array of time points to be saved while solving ODE. Typically defined as `t_save = t_span[1]:Δt_save:t_span[2]`
    - `dim_x`: `length(x)`
    - `dim_x_c`: `length(x_c)`

The keyword arguments should be provided as explained below:

- `solve_alg`: The algorithm used for solving ODEs. Default value is `Tsit5()`
- `sense_alg`: The algorithm used for adjoint sensitivity analysis. Default value is `InterpolatingAdjoint(autojacvec = ZygoteVJP())`, because the control problems usually render the `BacksolveAdjoint()` unstable. The vjp choice `autojacvec = ReverseDiffVJP(true)` is usually faster than `ZygoteVJP()`, when the ODE function does not have any branching inside. Please refer to the [DiffEqFlux documentation](https://diffeqflux.sciml.ai/dev/ControllingAdjoints/) for further details. 
- `ensemble_alg`: The algorithm used for handling ensemble of ODEs. Default value is `EnsembleThreads()` for multi-threaded computation in CPU.
- `opt_1`: The algorithm used for the first phase of optimisation which rapidly delivers the parameter to a favourable region around a local minimum. Default value is `ADAM(0.01)`.
- `opt_2`: The algorithm used for the second phase of opitmisaiton. Defalut value is `LBFGS()` which refines the result of the first phase to find a more precise minimum. Please refer to the [DiffEqFlux documentation](https://diffeqflux.sciml.ai/dev/sciml_train/) for further details about two-phase composition of optimisers.
- `maxiters_1`: The maximum number of iterations allowed for the first phase of optimisation with `opt_1`. Defalut value is `100`.
- `maxiters_2`: The maximum number of iterations allowed for the second phase of optimisation with `opt_2`. Defalut value is `100`.
- `progress_plot`: The indicator to plot the state history for a nominal condition among the ensemble during the learning process. Default value is `true`.
- `i_nominal`: The index to select the case to plot using `progress_plot` during optimisation process from the `ensemble` defined in `scenario`. Defalut value is `nothing`.
- `p_NN_0`: Initial value of the NN parameters supplied by the user to bypass random initialisation of `p_NN` or to continue optimisation from the previous result. Defalut value is `nothing`.
- `solve_kwargs...`: Additional keyword arguments that are passed onto the ODE solver.

`CTPG_train()` returns the following outputs:

- `result`: The final result of parameter optimisation.
- `fwd_ensemble_sol`: The ensemble solution of forward simulation using the final neural network parameters.
- `loss_history`: The history of loss function evaluated at each iteration.
"""

function CTPG_train(dynamics_plant::Function, dynamics_controller::Function, dynamics_sensor::Function, cost_running::Function, cost_terminal::Function, cost_regularisor::Function, policy_NN, scenario; solve_alg = Tsit5(), sense_alg = InterpolatingAdjoint(autojacvec = ZygoteVJP()), ensemble_alg = EnsembleThreads(), opt_1 = ADAM(0.01), opt_2 = LBFGS(), maxiters_1 = 100, maxiters_2 = 100, progress_plot = true, i_nominal = nothing, p_NN_0 = nothing, solve_kwargs...)

    # scenario parameters
    @unpack ensemble, t_span, t_save, dim_x, dim_x_c = scenario
    dim_ensemble = length(ensemble)
    mean_factor  = Float32(1 / dim_ensemble)
    dim_t_save   = length(t_save)
    if isnothing(i_nominal)
        i_nominal = max(round(Int, dim_ensemble / 2), 1)
    end

    # NN parameters initialisation
    if isnothing(p_NN_0)
        p_NN_0 = initial_params(policy_NN)
    end

    # augmented dynamics
    function fwd_dynamics(r)
        return function (x_aug, p_NN, t)
            x   = x_aug[1:dim_x]
            x_c = x_aug[dim_x+1:end-1]
            # ∫cost_running = x_aug[end]
            
            y            = dynamics_sensor(t, x)
            (dx_c, u, _) = dynamics_controller(t, x_c, y, r, p_NN, policy_NN)
            dx           = dynamics_plant(t, x, u)

            return [dx; dx_c; cost_running(t, x, y, u, r)]
        end
    end

    # ODE problem construction
    prob_base = ODEProblem(fwd_dynamics(ensemble[i_nominal].r), [ensemble[i_nominal].x₀; zeros(Float32, dim_x_c + 1)], t_span, p_NN_0)

    function generate_probs(p_NN)
        return function (prob, i, repeat)
            remake(prob, f = fwd_dynamics(ensemble[i].r), u0 = [ensemble[i].x₀; zeros(Float32, dim_x_c + 1)], p = p_NN)
        end
    end

    if dim_ensemble == 1
        prob_mtk = modelingtoolkitize(prob_base)
        prob_base = ODEProblem(prob_mtk, [], t_span, jac = true)
        ensemble_alg = EnsembleSerial()
    end

    # loss function definition
    function loss(p_NN)
        ensemble_prob = EnsembleProblem(prob_base, prob_func = generate_probs(p_NN))

        fwd_ensemble_sol_full = solve(ensemble_prob, solve_alg, ensemble_alg, saveat = t_save, trajectories = dim_ensemble, sensealg = sense_alg; solve_kwargs...)

        # version 1: mean(sol[end])
        # sol_length = Float32.([max(1,length(fwd_ensemble_sol_full[i])) for i in 1:dim_ensemble])
        # sol_length_ratio = maximum(sol_length)./sol_length

        # loss_val = mean([fwd_ensemble_sol_full[i][end][end] + cost_terminal(fwd_ensemble_sol_full[i][end][1:dim_x], ensemble[i].r) for i in 1:dim_ensemble] .* sol_length_ratio) + cost_regularisor(p_NN)
        
        # version 2: mean(Array[end])
        # loss_val = mean([(Array(fwd_ensemble_sol_full[i])[end,end] + cost_terminal(Array(fwd_ensemble_sol_full[i])[1:dim_x,end], ensemble[i].r)) for i in 1:dim_ensemble] .* sol_length_ratio) + cost_regularisor(p_NN)

        # version 3: scalar operation (best for code robustness as it can handle divergent case in ensemble)
        loss_val = 0.0f0
        for i in 1:dim_ensemble
            fwd_sol = Array(fwd_ensemble_sol_full[i])

            if size(fwd_sol,2) > dim_t_save
                fwd_sol = fwd_sol[:,1:dim_t_save]
            end

            if size(fwd_sol,2) == dim_t_save
                x_aug_f = fwd_sol[:,end]
                loss_val += (x_aug_f[end] + cost_terminal(x_aug_f[1:dim_x], ensemble[i].r)) * mean_factor
            else
                loss_val += 1000.0f0
            end
        end
        loss_val += cost_regularisor(p_NN)
        
        return loss_val, fwd_ensemble_sol_full
    end

    # learning progress callback setup
    loss_history = fill(NaN32, maxiters_1 + maxiters_2 + 2)
    iterator_learning = 1
    cb_progress = function (p_NN_val, loss_val, fwd_ensemble_sol_full; plot_val = progress_plot)
        @show (loss_val, iterator_learning);
        loss_history[iterator_learning] = loss_val
        if plot_val
            fwd_ensemble_sol = Array(fwd_ensemble_sol_full[i_nominal])
            display(scatter(fwd_ensemble_sol[1:dim_x, :]', label = :false, plot_title = "System State: Learning Iteration $(iterator_learning)", layout = (dim_x, 1), size = (700, 200 * dim_x)))
        end
        iterator_learning += 1
        return false
    end

    # NN training
    if maxiters_1 <= 0
        result_coarse = (; u = p_NN_0)
    else
        result_coarse = DiffEqFlux.sciml_train(loss, p_NN_0, opt_1; cb = cb_progress, maxiters = maxiters_1)
    end
    
    if maxiters_2 <= 0
        result = result_coarse
    else
        result = DiffEqFlux.sciml_train(loss, result_coarse.u, opt_2; cb = cb_progress, maxiters = maxiters_2)
    end

    # Forward solution for optimised p_NN 
    function eval_IO(r, p_NN)
        # Evaluation of input u, output y, and NN output y_NN at each t
        return function (t, x_aug)
            x = x_aug[1:dim_x]
            x_c = x_aug[dim_x+1:end-1]

            y = dynamics_sensor(t, x)
            (_, u, y_NN) = dynamics_controller(t, x_c, y, r, p_NN, policy_NN)

            return (; u = u, y = y, y_NN = y_NN)
        end
    end

    function generate_savedata(p_NN)
        return function (sol, i)
            IO = eval_IO(ensemble[i].r, p_NN).(sol.t, sol.u)
            return ((; sol = sol, u = [u for (u, y, y_NN) in IO], y = [y for (u, y, y_NN) in IO], y_NN = [y_NN for (u, y, y_NN) in IO]), false)
        end
    end

    fwd_ensemble_sol = solve(EnsembleProblem(prob_base, prob_func = generate_probs(result.u), output_func = generate_savedata(result.u)), solve_alg, ensemble_alg, saveat = t_save, trajectories = dim_ensemble, sensealg = sense_alg; solve_kwargs...)

    loss_history = filter(!isnan, loss_history)

    return result, fwd_ensemble_sol, loss_history
end
