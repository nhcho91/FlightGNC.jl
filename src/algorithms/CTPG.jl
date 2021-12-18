function CTPG_train(dynamics_plant::Function, dynamics_controller::Function, cost_running::Function, cost_terminal::Function, cost_regularisor::Function, policy_NN, scenario; 
    solve_alg = Tsit5(), sense_alg = InterpolatingAdjoint(autojacvec = ZygoteVJP()), ensemble_alg = EnsembleThreads(), opt_1 = ADAM(0.01), opt_2 = LBFGS(), progress_plot = true, solve_kwargs...)
    
    # scenario parameters
    @unpack ensemble_list, t_span, dim_x, dim_x_c = scenario
    dim_ensemble = length(ensemble_list)
    id_nominal = max(round(Int, dim_ensemble / 2), 1)

    # NN parameters initialisation
    p_NN_0 = initial_params(policy_NN)

    # augmented dynamics
    function fwd_dynamics(r)
        return function (x_aug, p_NN, t)
            x   = x_aug[1:dim_x]
            x_c = x_aug[dim_x+1:end-1]
            # ∫cost_running = x_aug[end]

            (dx_c, u) = dynamics_controller(t, x_c, x, r, p_NN, policy_NN)
            (dx, y) = dynamics_plant(t, x, u)

            return [dx; dx_c; cost_running(t, x, y, u, r)]
        end
    end

    # ODE problem construction
    prob_base = ODEProblem(fwd_dynamics(ensemble_list[id_nominal][2]), [ensemble_list[id_nominal][1]; zeros(Float32, dim_x_c + 1)], t_span, p_NN_0)

    function generate_probs(p_NN)
        return function (prob, i, repeat)
            remake(prob, f = fwd_dynamics(ensemble_list[i][2]), u0 = [ensemble_list[i][1]; zeros(Float32, dim_x_c + 1)], p = p_NN)
        end
    end

    # loss function definition
    function loss(p_NN)
        ensemble_prob = EnsembleProblem(prob_base, prob_func = generate_probs(p_NN))

        fwd_ensemble_sol = Array(solve(ensemble_prob, solve_alg, ensemble_alg, trajectories = dim_ensemble, sensealg = sense_alg; solve_kwargs...))

        loss_val = mean(cost_terminal(fwd_ensemble_sol[1:dim_x, end, i], ensemble_list[i][2]) for i = 1:dim_ensemble) + mean(fwd_ensemble_sol[end, end, :]) + cost_regularisor(p_NN)

        fwd_sol_nominal_val = fwd_ensemble_sol[:, :, id_nominal]

        return loss_val, fwd_sol_nominal_val
    end

    # learning progress callback setup
    loss_history = fill(NaN32, 250)
    iterator_learning = 1
    cb_progress = function (p_NN_val, loss_val, fwd_sol_nominal_val; plot_val = progress_plot)
        @show loss_val
        loss_history[iterator_learning] = loss_val
        if plot_val
            display(scatter(fwd_sol_nominal_val[1:dim_x, :]', label = :false, plot_title = "System State: Learning Iteration $(iterator_learning)", layout = (dim_x, 1), size = (700, 200 * dim_x)))
        end
        iterator_learning += 1
        return false
    end

    # NN training
    result_coarse = DiffEqFlux.sciml_train(loss, p_NN_0,          opt_1; cb = cb_progress, maxiters = 50)
    result        = DiffEqFlux.sciml_train(loss, result_coarse.u, opt_2; cb = cb_progress, maxiters = 50)
    # result = DiffEqFlux.sciml_train(loss, p_NN_0; cb = cb_progress, maxiters = 100)

    fwd_sol_nominal = solve(prob_base, solve_alg, u0 = [ensemble_list[id_nominal][1]; zeros(Float32, dim_x_c + 1)], p = result.u; solve_kwargs...)

    return result, fwd_sol_nominal, loss_history
end
