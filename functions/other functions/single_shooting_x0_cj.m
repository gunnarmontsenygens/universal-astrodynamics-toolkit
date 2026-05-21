function z_vec = single_shooting_x0_cj(z_0_vec, g_i_hat, Cj_star, params, iter_max, tol, show)

    % Initialization
    z_vec = z_0_vec(:);
    g_i_hat = g_i_hat(:).';

    % Build selection matrix S by removing the row selected by g_i_hat
    idx = find(abs(g_i_hat) > 0, 1);

    if isempty(idx)
        error('g_i_hat must select one state component.');
    end

    S = eye(6);
    S(idx,:) = [];

    % Model selection
    switch upper(params.model.name)
        case 'CR3BP'
            eom = @eom_cr3bp;
            integrate = @integrate_cr3bp;
            jacobi_constant = @jacobi_constant_cr3bp;
            jacobi_gradient = @jacobi_gradient_cr3bp;

        case 'HILLR3BP'
            eom = @eom_hillr3bp;
            integrate = @integrate_hillr3bp;
            jacobi_constant = @jacobi_constant_hillr3bp;
            jacobi_gradient = @jacobi_gradient_hillr3bp;

        otherwise
            error('Model not registered')
    end

    for i = 1:iter_max

        % -------------------------------------------------------------
        % Step I: EXTRACT x0_vec AND T
        % -------------------------------------------------------------
        z_0_vec = z_vec;
        x_0_vec = z_0_vec(1:6);
        T = z_0_vec(7);

        % -------------------------------------------------------------
        % Step II: INTEGRATE FOR T TIME UNITS
        % -------------------------------------------------------------
        [~, x_vec_hist, Phi_mtx_hist, ~, ~, ~, ~] = ...
            integrate([0,T], x_0_vec, params);

        % -------------------------------------------------------------
        % STEP III: CORRECT
        % -------------------------------------------------------------

        % Extract final state
        x_T_vec = x_vec_hist(end,:)';

        % Extract final STM
        Phi_T_mtx = squeeze(Phi_mtx_hist(end,:,:));

        % Vector field at final state
        dx_dt_T_vec = eom(T, x_T_vec, params);
        dx_dt_T_vec = dx_dt_T_vec(1:6);
        dx_dt_T_vec = dx_dt_T_vec(:);

        % Jacobi gradient at initial state
        dCjdx_0_vec = jacobi_gradient(x_0_vec, params).';

        % Constraint vector
        g = g_i_hat*x_0_vec;

        Cj_0 = jacobi_constant(0, x_0_vec.', params);

        h_vec = [S*(x_T_vec - x_0_vec);
                 g;
                 Cj_0 - Cj_star];

        % Jacobian wrt z = [x0; T]
        dhdz_mtx = [S*(Phi_T_mtx - eye(6)), S*dx_dt_T_vec;
                    g_i_hat,                0;
                    dCjdx_0_vec,            0];

        % Newton correction
        Delta_z_vec = -pinv(dhdz_mtx)*h_vec;

        % Maximum allowed relative period correction
        dT_max = 0.2*abs(T);

        % Limit period correction magnitude
        if abs(Delta_z_vec(7)) > dT_max

            % Scale period correction
            Delta_z_vec(7) = sign(Delta_z_vec(7))*dT_max;

        end

        % -------------------------------------------------------------
        % STEP V: DAMPED NEWTON CORRECTION, LINE SEARCH AND UPDATE
        % -------------------------------------------------------------

        % Initialize full Newton step
        alpha = 1.0;
        
        % Minimum allowable damping factor
        alpha_min = 1e-6;
        
        % Store current constraint norm
        h_norm_old = norm(h_vec);
        
        % Begin backtracking line search
        while alpha > alpha_min
        
            % Compute trial correction step
            z_trial_vec = z_0_vec + alpha*Delta_z_vec;
        
            % Minimum allowable period
            T_min = 0.5;

            % Maximum allowable period
            T_max = 20.0;

            % Check if trial period is outside allowable bounds
            if z_trial_vec(7) < T_min || z_trial_vec(7) > T_max
        
                % Reduce step size
                alpha = 0.5*alpha;
        
                % Skip remaining loop body
                continue
        
            end
        
            % Extract corrected initial state
            x_trial_0_vec = z_trial_vec(1:6);
        
            % Extract corrected period
            T_trial = z_trial_vec(7);
        
            % Integrate trajectory using trial correction
            [~, x_trial_hist, ~, ~, ~, ~, ~] = ...
                integrate([0,T_trial], x_trial_0_vec, params);
        
            % Extract final state after propagation
            x_trial_T_vec = x_trial_hist(end,:)';
        
            % Compute Jacobi constant at corrected initial condition
            Cj_trial = jacobi_constant(0, x_trial_0_vec.', params);
        
            % Compute trial constraint vector
            h_trial_vec = [S*(x_trial_T_vec - x_trial_0_vec);
                           g_i_hat*x_trial_0_vec;
                           Cj_trial - Cj_star];
        
            % Check whether trial step improves the solution
            if norm(h_trial_vec) < h_norm_old
        
                % Accept current damping factor
                break
        
            end
        
            % Reduce Newton step size
            alpha = 0.5*alpha;
        
        end

        % Check if line search failed
        if alpha <= alpha_min
        
            % Stop correction if no acceptable step was found
            warning('Line search failed at iteration %d.', i)
        
            % Exit Newton loop
            break
        
        end
        
        % Update solution using accepted damped Newton step
        z_vec = z_0_vec + alpha*Delta_z_vec;

        % -------------------------------------------------------------
        % STEP IV: CHECK TOLERANCE
        % -------------------------------------------------------------
        h = norm(h_vec);
        
        if show
            fprintf('Iteration %d. h = %.3e, dz = %.3e\n', ...
            i, h, norm(Delta_z_vec))
        end

        if h < tol
            break
        end

    end

end