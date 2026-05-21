function z_vec = single_shooting_x0(z_0_vec, g_i_hat, params, iter_max, tol, show)
%==========================================================================
%
% Computes a periodic orbit using a single-shooting differential correction
% scheme for a specified dynamical model within the Universal Modeling
% Toolkit framework.
%
% METHOD DESCRIPTION:
% The algorithm solves for an initial condition and period satisfying the
% periodicity condition:
%
%       x(T) = x_0
%
% together with an initial-state phase/symmetry constraint:
%
%       g_i_hat * x_0 = 0
%
% The correction problem is solved using a damped Newton iteration with a
% backtracking line search. The State Transition Matrix (STM) is used to
% construct the Jacobian of the shooting constraints.
%
% UNKNOWN VECTOR:
%
%       z = [x_0; T]
%
% where:
%   x_0     initial state vector
%   T       orbit period
%
% CONSTRAINT VECTOR:
%
%       h(z) =
%       [x(T;x_0) - x_0;
%        g_i_hat*x_0]
%
% The system is solved iteratively using a minimum-norm Newton correction.
%
% NOTES:
% - The initial-state phase condition removes the time-shift degeneracy
%   of periodic orbits.
% - A damped Newton correction and line search are used to improve
%   robustness near unstable or highly nonlinear solutions.
% - The algorithm supports multiple dynamical models through the model
%   definition stored in params.model.
%
% Author: G. Montseny
% Date: May 7, 2026
%
% INPUT:               Description                                   Units
%
%  z_0_vec   -   initial guess [x_0; T]                              [-]
%  g_i_hat   -   phase/symmetry selection row vector                 [-]
%  params    -   parameter struct defining the model                 [-]
%  iter_max  -   maximum number of Newton iterations                 [-]
%  tol       -   convergence tolerance                               [-]
%  show      -   logical flag for iteration output                   [-]
%
% OUTPUT:              Description                                   Units
%
%  z_vec     -   corrected periodic orbit state and period           [-]
%
%==========================================================================

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

        case 'HILLR3BP'
            eom = @eom_hillr3bp;
            integrate = @integrate_hillr3bp;

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


        % Constraint vector
        g = g_i_hat*x_0_vec;

        h_vec = [(x_T_vec - x_0_vec);
                 g];

        % Jacobian wrt z = [x0; T]
        dhdz_mtx = [(Phi_T_mtx - eye(6)), dx_dt_T_vec;
                    g_i_hat,                0];

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
        
        
            % Compute trial constraint vector
            h_trial_vec = [S*(x_trial_T_vec - x_trial_0_vec);
                           g_i_hat*x_trial_0_vec];
        
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