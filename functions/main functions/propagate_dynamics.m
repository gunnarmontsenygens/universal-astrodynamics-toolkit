function traj = propagate_dynamics(t_span, x_0_vec, params, event_fun, symp_tol)
%==========================================================================
%
% Propagates the state and State Transition Matrix (STM) for a specified
% dynamical model within the Universal Modeling Toolkit framework. This
% function serves as a high-level interface that dispatches the integration
% to the appropriate model-specific integrator based on the model definition
% contained in the parameter structure.
%
% MODEL DESCRIPTION:
% The dynamical system is defined through the fields in params.model, which
% specify:
%   - model.name         (e.g., 'CR3BP', 'CCR4BP', 'NBP')
%   - model.formulation  (e.g., 'lagrangian', 'hamiltonian')
%   - model.frame        (e.g., 'synodic', 'inertial')
%   - model.units        (e.g., 'nd', 'dim')
%
% Based on this configuration, the appropriate low-level integrator is
% called (e.g., integrate_cr3bp for the CR3BP in Lagrangian synodic
% nondimensional form).
%
% STATE DEFINITION:
% The propagated state corresponds to the augmented system:
%
%   X_vec = [x; vec(Phi)]
%
% where:
%   - x is the physical state vector
%   - Phi is the State Transition Matrix (STM)
%
% OUTPUT STRUCTURE:
% The output is a trajectory structure 'traj' containing:
%
%   traj.model          - model definition (params.model)
%   traj.params         - full parameter struct
%   traj.t_span         - integration interval
%   traj.x_0_vec        - initial state
%   traj.t_hist         - time history
%   traj.x_vec_hist     - state history
%   traj.Phi_mtx_hist   - STM history
%   traj.x_f_vec        - final state
%   traj.Phi_f_mtx      - final STM
%
%   traj.events         - event data:
%       i_e             - event indices
%       t_e             - event times
%       x_vec           - state at events
%       Phi_mtx         - STM at events
%
%   traj.invariants     - model-dependent invariant quantities
%       (e.g., energy, Jacobi constant for CR3BP)
%
%   traj.symp           - symplecticity diagnostics 
%       t_hist          - time history used for symplectic check
%       err_hist        - symplectic error history (||Phi^T J Phi - J||_F)
%       err_f           - final symplectic error
%       err_max         - maximum symplectic error over trajectory
%       valid           - logical flag (true if err_max < tol)
%       tol             - tolerance used for symplecticity check
%       formulation_checked - formulation used for the check (if converted)
%       converted_from  - original formulation if conversion was applied
%
% NOTES:
% - This function separates the dynamical model definition (params) from the
%   numerical experiment (event_fun), ensuring modularity and extensibility.
% - Additional models (e.g., CCR4BP, NBP) can be incorporated by extending
%   the dispatch structure.
% - If the system is not provided in Hamiltonian form, it is internally
%   converted before performing the symplecticity check.
%
%
% Author: G. Montseny
% Date: May 4, 2026
%
% INPUT:               Description                                   Units
%
%  t_span     -   time span for integration [t0 tf]                  [-]
%  x_0_vec    -   initial state vector                               [-]
%  params     -   parameter struct defining model and constants      [-]
%  event_fun  -   optional event function handle                     [-]
%  symp_tol   -   optional tolerance for symplecticity check         [-]
%
% OUTPUT:              Description                                   Units
%
%  traj       -   trajectory structure containing state, STM,        [-]
%                 events, invariant quantities, and symplectic data
%
%==========================================================================

    % Default event function
    if nargin < 4
        event_fun = [];
    end

    % Default symplecticity tolerance

    if nargin < 5 || isempty(symp_tol)
        symp_tol = 1e-10;
    end

    % Extract model identity
    model_name  = upper(params.model.name);
    formulation = lower(params.model.formulation);
    frame       = lower(params.model.frame);
    units       = lower(params.model.units);

    % Dispatch to model-specific integrator
    switch model_name

        case 'CR3BP'

            if strcmp(formulation,'lagrangian') && strcmp(frame,'synodic') && strcmp(units,'nd')
                [t_hist, x_vec_hist, Phi_mtx_hist, i_e, t_e, x_e_vec, Phi_mtx_e] = integrate_cr3bp(t_span, x_0_vec, params, event_fun);
            else
                error('Unsupported CR3BP case: %s, %s, %s.', formulation, frame, units);
            end

        case 'HILLR3BP'

            if strcmp(formulation,'lagrangian') && strcmp(frame,'synodic') && strcmp(units,'nd')
                [t_hist, x_vec_hist, Phi_mtx_hist, i_e, t_e, x_e_vec, Phi_mtx_e] = integrate_hillr3bp(t_span, x_0_vec, params, event_fun);
            else
                error('Unsupported HILLR3BP case: %s, %s, %s.', formulation, frame, units);
            end

        otherwise
            error('Unsupported model: %s.', model_name);

    end

    % Build trajectory struct
    traj = struct();

    traj.model = params.model;
    traj.params = params;

    traj.t_span = t_span;
    traj.x_0_vec = x_0_vec(:);

    traj.t_hist = t_hist;
    traj.x_vec_hist = x_vec_hist;
    traj.Phi_mtx_hist = Phi_mtx_hist;

    traj.x_f_vec = x_vec_hist(end,:)';
    traj.Phi_f_mtx = squeeze(Phi_mtx_hist(end,:,:));

    % Events
    traj.events.i_e = i_e;
    traj.events.t_e = t_e;
    traj.events.x_vec = x_e_vec;
    traj.events.Phi_mtx = Phi_mtx_e;

    % Invariants and symplecticity 
    switch model_name
        case 'CR3BP'

            traj.invariants.E_hist = energy_cr3bp(t_hist, x_vec_hist, params);
            traj.invariants.C_hist = jacobi_constant_cr3bp(t_hist, x_vec_hist, params);
            traj.invariants.DeltaC_hist = traj.invariants.C_hist - traj.invariants.C_hist(1);
    
            if strcmp(formulation,'hamiltonian')
                [traj.symp.t_hist, ...
                 traj.symp.err_hist, ...
                 traj.symp.err_f, ...
                 traj.symp.err_max, ...
                 traj.symp.valid] = symp_cr3bp(t_hist, Phi_mtx_hist, symp_tol, params);
    
            else
    
                [t_hist_h, ~, Phi_mtx_hist_h, ~, ~, ~, ~, params_h] = ...
                    l2h_cr3bp(t_hist, x_vec_hist, Phi_mtx_hist, i_e, t_e, x_e_vec, Phi_mtx_e, params);
                [traj.symp.t_hist, ...
                 traj.symp.err_hist, ...
                 traj.symp.err_f, ...
                 traj.symp.err_max, ...
                 traj.symp.valid] = symp_cr3bp(t_hist_h, Phi_mtx_hist_h, symp_tol, params_h);
    
                traj.symp.formulation_checked = 'hamiltonian';
                traj.symp.converted_from = formulation;
    
            end
    
            traj.symp.tol = symp_tol;

        case 'HILLR3BP'

            traj.invariants.J_hist = jacobi_integral_hillr3bp(t_hist, x_vec_hist, params);
            traj.invariants.C_hist = jacobi_constant_hillr3bp(t_hist, x_vec_hist, params);
            traj.invariants.DeltaC_hist = traj.invariants.C_hist - traj.invariants.C_hist(1);
    
            if strcmp(formulation,'hamiltonian')
                [traj.symp.t_hist, ...
                 traj.symp.err_hist, ...
                 traj.symp.err_f, ...
                 traj.symp.err_max, ...
                 traj.symp.valid] = symp_hillr3bp(t_hist, Phi_mtx_hist, symp_tol, params);
    
            else
    
                [t_hist_h, ~, Phi_mtx_hist_h, ~, ~, ~, ~, params_h] = ...
                    l2h_hillr3bp(t_hist, x_vec_hist, Phi_mtx_hist, i_e, t_e, x_e_vec, Phi_mtx_e, params);
                [traj.symp.t_hist, ...
                 traj.symp.err_hist, ...
                 traj.symp.err_f, ...
                 traj.symp.err_max, ...
                 traj.symp.valid] = symp_hillr3bp(t_hist_h, Phi_mtx_hist_h, symp_tol, params_h);
    
                traj.symp.formulation_checked = 'hamiltonian';
                traj.symp.converted_from = formulation;
    
            end
    
            traj.symp.tol = symp_tol;

        otherwise
            error('Model name not recognized.')
    end
end