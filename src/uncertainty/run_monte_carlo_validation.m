function mc_results = run_monte_carlo_validation(sim_params, R_intensity, duration_s, pse_config, varargin)
% RUN_MONTE_CARLO_VALIDATION  Monte Carlo validation of PSE uncertainty bounds
%
% Validates DAE-based PSE uncertainty bounds using independent solvers
% (CA_Routing or semi-implicit). Samples parameters from the same
% distributions PSE uses, runs N realizations, and compares ensemble
% statistics against PSE predictions.
%
% KEY DESIGN: PSE propagates uncertainty through the DAE Jacobian.
% MC validation uses INDEPENDENT solvers (explicit CA or semi-implicit)
% to confirm PSE bounds are physically valid, not just DAE-consistent.
%
% INPUTS:
%   sim_params   - Domain parameters (from setup functions)
%   R_intensity  - Rainfall forcing (struct or scalar)
%   duration_s   - Simulation duration [s]
%   pse_config   - PSE configuration with CV values
%
% NAME-VALUE PAIRS:
%   'Solver'          - 'ca_routing' (default) or 'semi_implicit'
%   'NSamples'        - Number of MC samples (default: 100)
%   'Sampling'        - 'lhs' (default) or 'random'
%   'UseParallel'     - Use parfor if available (default: true)
%   'SaveTrajectories'- Save all N trajectories (default: false)
%   'DisplayProgress' - Print progress (default: true)
%   'SolverDt'        - Solver time step [s] (default: 1.0 for CA, 0.5 for SI)
%   'OutputInterval'  - CA output interval [steps] (default: 60)
%   'RandomSeed'      - RNG seed for reproducibility (default: 42)
%
% OUTPUTS:
%   mc_results - Structure with:
%     .h_mean_final      - Mean final depth [ny x nx]
%     .h_std_final        - Std final depth [ny x nx]
%     .h_p025_final       - 2.5th percentile final depth [ny x nx]
%     .h_p975_final       - 97.5th percentile final depth [ny x nx]
%     .Q_outlet_ensemble  - Outlet Q time series [N x nt] (if SaveTrajectories)
%     .Q_outlet_mean      - Mean outlet Q [nt x 1]
%     .Q_outlet_std       - Std outlet Q [nt x 1]
%     .Q_outlet_p025      - 2.5th percentile Q [nt x 1]
%     .Q_outlet_p975      - 97.5th percentile Q [nt x 1]
%     .t_common           - Common time vector [nt x 1]
%     .F_mean_final       - Mean final infiltration [ny x nx]
%     .F_std_final         - Std final infiltration [ny x nx]
%     .param_samples       - Sampled parameter values [N x n_params]
%     .param_names         - Cell array of parameter names
%     .n_samples           - Number of samples
%     .n_successful        - Number of successful runs
%     .solver_name         - Solver used
%     .random_seed         - Random seed for reproducibility
%     .sampling_method     - Sampling method ('lhs' or 'random')
%     .computation_time    - Total wall clock [s]
%     .convergence         - Running σ estimate at 10%, 25%, 50%, 75%, 100%
%
% REFERENCE:
%   Validates  PSE bounds using solver-independent MC
%
% Author: M. H. Kazma
% Date: February 2026

%% Parse inputs
p = inputParser;
addRequired(p, 'sim_params', @isstruct);
addRequired(p, 'R_intensity');
addRequired(p, 'duration_s', @isnumeric);
addRequired(p, 'pse_config', @isstruct);
addParameter(p, 'Solver', 'ca_routing', @(x) ismember(x, {'ca_routing', 'semi_implicit'}));
addParameter(p, 'NSamples', 100, @isnumeric);
addParameter(p, 'Sampling', 'lhs', @(x) ismember(x, {'lhs', 'random'}));
addParameter(p, 'UseParallel', true, @islogical);
addParameter(p, 'SaveTrajectories', false, @islogical);
addParameter(p, 'DisplayProgress', true, @islogical);
addParameter(p, 'SolverDt', [], @isnumeric);
addParameter(p, 'OutputInterval', 60, @isnumeric);
addParameter(p, 'RandomSeed', 42, @isnumeric);
addParameter(p, 'LightweightSemi', false, @islogical);  % Skip 3D history in semi-implicit (saves RAM)
parse(p, sim_params, R_intensity, duration_s, pse_config, varargin{:});

solver_name     = p.Results.Solver;
N               = p.Results.NSamples;
sampling_method = p.Results.Sampling;
use_parallel    = p.Results.UseParallel;
save_traj       = p.Results.SaveTrajectories;
display_prog    = p.Results.DisplayProgress;
output_interval = p.Results.OutputInterval;
random_seed     = p.Results.RandomSeed;
lightweight_semi = p.Results.LightweightSemi;

% Set random seed for reproducibility
if ~isempty(random_seed) && isnumeric(random_seed)
    rng(random_seed, 'twister');
    if display_prog
        fprintf('  Random seed: %d (reproducible)\n', random_seed);
    end
end

% Default dt per solver
if isempty(p.Results.SolverDt)
    if strcmp(solver_name, 'ca_routing')
        solver_dt = 1.0;
    else
        solver_dt = 0.5;
    end
else
    solver_dt = p.Results.SolverDt;
end

ny = sim_params.grid_size(1);
nx = sim_params.grid_size(2);

fprintf('\n=== Monte Carlo Validation ===\n');
fprintf('  Solver: %s\n', solver_name);
fprintf('  Samples: %d (%s)\n', N, sampling_method);
fprintf('  Domain: %d x %d cells\n', ny, nx);
fprintf('  Solver dt: %.2f s\n', solver_dt);

tic_total = tic;

%% ========================================================================
%  STEP 1: Identify active parameters and sample
% =========================================================================
[param_names, param_base, param_sigma] = identify_active_params(sim_params, pse_config);
n_params = numel(param_names);

fprintf('  Active uncertainty sources: %d\n', n_params);
for pp = 1:n_params
    fprintf('    %s: base=%.4e, sigma=%.4e (CV=%.0f%%)\n', ...
        param_names{pp}, param_base(pp), param_sigma(pp), ...
        100 * param_sigma(pp) / max(param_base(pp), eps));
end

% Generate samples
if strcmp(sampling_method, 'lhs')
    U = lhsdesign(N, n_params);
else
    U = rand(N, n_params);
end

% Get parameter distributions from PSE config (match PSE assumptions)
param_dist = get_param_distributions(pse_config, param_names);

% Map uniform samples to parameter distributions according to PSE config
param_samples = zeros(N, n_params);
for pp = 1:n_params
    mu = param_base(pp);
    sigma = param_sigma(pp);
    dist_type = param_dist{pp};
    
    if mu > 0 && sigma > 0
        if strcmpi(dist_type, 'lognormal')
            % Lognormal: match mean and variance to mu, sigma²
            mu_ln = log(mu^2 / sqrt(sigma^2 + mu^2));
            sigma_ln = sqrt(log(1 + sigma^2 / mu^2));
            param_samples(:, pp) = exp(mu_ln + sigma_ln * norminv(U(:, pp)));
        elseif strcmpi(dist_type, 'truncated_gaussian')
            % Truncated Gaussian: sample then clip to [0, 1]
            samples = mu + sigma * norminv(U(:, pp));
            param_samples(:, pp) = max(0, min(1, samples));
        else
            % Default: Gaussian
            param_samples(:, pp) = mu + sigma * norminv(U(:, pp));
        end
    else
        param_samples(:, pp) = mu * ones(N, 1);
    end
end

%% ========================================================================
%  STEP 2: Run ensemble
% =========================================================================
fprintf('  Running %d realizations...\n', N);

% Pre-allocate per-realization collectors
h_final_all = zeros(N, ny, nx);
F_final_all = zeros(N, ny, nx);
success_flags = false(N, 1);

% Common time grid for outlet Q (estimated from single run schedule)
% Will be refined after first successful run
t_common = [];
Q_outlet_all = [];

% Convergence checkpoints
check_fracs = [0.10, 0.25, 0.50, 0.75, 1.00];
convergence = struct('frac', check_fracs, 'h_std_mean', nan(size(check_fracs)));

% Check parallel availability
has_parallel = use_parallel && ~isempty(ver('parallel'));
if has_parallel
    try
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool;
        end
        fprintf('  Parallel pool: %d workers\n', pool.NumWorkers);
    catch
        has_parallel = false;
        fprintf('  Parallel pool unavailable, using serial.\n');
    end
end

if has_parallel
    % === PARFOR Execution ===
    % Pre-build perturbed sim_params for each sample
    parfor_results = cell(N, 1);

    parfor m = 1:N
        try
            % Perturb parameters
            sp_m = perturb_params(sim_params, param_names, param_samples(m, :)); %#ok<PFBNS>
            R_m = perturb_rainfall(R_intensity, pse_config, param_names, param_samples(m, :));

            % Run solver
            if strcmp(solver_name, 'ca_routing')
                res_m = run_ca_routing_comparison(sp_m, R_m, duration_s, ...
                    'dt', solver_dt, 'UseInfiltration', true, ...
                    'OutputInterval', output_interval, 'DisplayProgress', false);
                parfor_results{m} = struct('h_final', res_m.h_final, ...
                    'F_final', res_m.F_final, ...
                    'Q_outlet', res_m.Q_outlet, ...
                    't', res_m.t, 'ok', true);
            else
                [t_m, res_m] = run_semi_implicit_simulation(sp_m, R_m, duration_s, ...
                    'dt', solver_dt, 'DisplayProgress', false, ...
                    'LightweightMode', lightweight_semi);
                % In lightweight mode, res_m.h is (ny,nx) final only;
                % in normal mode it is (ny,nx,nt) history.
                if ndims(res_m.h) == 3 %#ok<ISMAT>
                    h_f = res_m.h(:,:,end);
                else
                    h_f = res_m.h;
                end
                F_f = zeros(ny, nx);
                if isfield(res_m, 'F')
                    if ndims(res_m.F) == 3 %#ok<ISMAT>
                        F_f = res_m.F(:,:,end);
                    else
                        F_f = res_m.F;
                    end
                end
                % FIX [AWR2026]: Use Q_outlet_ts (flume/designated cells
                % only) to match CA Routing Q_outlet measurement.
                Q_out = [];
                if isfield(res_m, 'Q_outlet_ts')
                    Q_out = res_m.Q_outlet_ts;
                elseif isfield(res_m, 'Q_boundary_ts')
                    Q_out = res_m.Q_boundary_ts;  % fallback
                end
                parfor_results{m} = struct('h_final', h_f, 'F_final', F_f, ...
                    'Q_outlet', Q_out, 't', t_m, 'ok', true);
            end
        catch
            parfor_results{m} = struct('h_final', [], 'F_final', [], ...
                'Q_outlet', [], 't', [], 'ok', false);
        end
    end

    % Unpack
    for m = 1:N
        if parfor_results{m}.ok
            success_flags(m) = true;
            h_final_all(m, :, :) = parfor_results{m}.h_final;
            F_final_all(m, :, :) = parfor_results{m}.F_final;
            if isempty(t_common) && ~isempty(parfor_results{m}.t)
                t_common = parfor_results{m}.t(:);
                Q_outlet_all = nan(N, numel(t_common));
            end
            if ~isempty(parfor_results{m}.Q_outlet) && ~isempty(t_common)
                Q_m = parfor_results{m}.Q_outlet(:);
                t_m = parfor_results{m}.t(:);
                if numel(Q_m) == numel(t_m)
                    Q_interp = interp1(t_m, Q_m, t_common, 'linear', 'extrap');
                    Q_outlet_all(m, :) = Q_interp;
                end
            end
        end
    end

    if display_prog
        fprintf('  Parallel MC complete: %d/%d successful\n', sum(success_flags), N);
    end

else
    % === SERIAL Execution ===
    for m = 1:N
        try
            sp_m = perturb_params(sim_params, param_names, param_samples(m, :));
            R_m = perturb_rainfall(R_intensity, pse_config, param_names, param_samples(m, :));

            if strcmp(solver_name, 'ca_routing')
                res_m = run_ca_routing_comparison(sp_m, R_m, duration_s, ...
                    'dt', solver_dt, 'UseInfiltration', true, ...
                    'OutputInterval', output_interval, 'DisplayProgress', false);
                h_f = res_m.h_final;
                F_f = res_m.F_final;
                Q_out = res_m.Q_outlet;
                t_m = res_m.t;
            else
                [t_m, res_m] = run_semi_implicit_simulation(sp_m, R_m, duration_s, ...
                    'dt', solver_dt, 'DisplayProgress', false, ...
                    'LightweightMode', lightweight_semi);
                if ndims(res_m.h) == 3 %#ok<ISMAT>
                    h_f = res_m.h(:,:,end);
                else
                    h_f = res_m.h;
                end
                F_f = zeros(ny, nx);
                if isfield(res_m, 'F')
                    if ndims(res_m.F) == 3 %#ok<ISMAT>
                        F_f = res_m.F(:,:,end);
                    else
                        F_f = res_m.F;
                    end
                end
                % FIX [AWR2026]: Use Q_outlet_ts (flume/designated cells
                % only) to match CA Routing Q_outlet measurement.
                Q_out = [];
                if isfield(res_m, 'Q_outlet_ts')
                    Q_out = res_m.Q_outlet_ts;
                elseif isfield(res_m, 'Q_boundary_ts')
                    Q_out = res_m.Q_boundary_ts;  % fallback
                end
            end

            success_flags(m) = true;
            h_final_all(m, :, :) = h_f;
            F_final_all(m, :, :) = F_f;

            % Initialize common time grid from first success
            if isempty(t_common) && ~isempty(t_m)
                t_common = t_m(:);
                Q_outlet_all = nan(N, numel(t_common));
            end
            if ~isempty(Q_out) && ~isempty(t_common)
                Q_m = Q_out(:);
                if numel(Q_m) == numel(t_m(:))
                    Q_interp = interp1(t_m(:), Q_m, t_common, 'linear', 'extrap');
                    Q_outlet_all(m, :) = Q_interp;
                end
            end

        catch ME_m
            if display_prog
                fprintf('  Sample %d FAILED: %s\n', m, ME_m.message);
            end
        end

        % Progress and convergence check
        if display_prog && (mod(m, max(1, round(N/10))) == 0 || m == N)
            n_ok = sum(success_flags(1:m));
            if n_ok > 1
                h_ok = h_final_all(success_flags(1:m), :, :);
                running_std = mean(std(h_ok, 0, 1), 'all');
                fprintf('  [%3d/%d] ok=%d, running mean σ_h=%.4f mm\n', ...
                    m, N, n_ok, running_std * 1000);

                % Store convergence
                for ci = 1:numel(check_fracs)
                    if m >= round(check_fracs(ci) * N) && isnan(convergence.h_std_mean(ci))
                        convergence.h_std_mean(ci) = running_std;
                    end
                end
            end
        end
    end
end

%% ========================================================================
%  STEP 3: Compute ensemble statistics
% =========================================================================
n_ok = sum(success_flags);
fprintf('  Computing statistics from %d successful runs...\n', n_ok);

if n_ok < 2
    error('Fewer than 2 successful MC runs. Cannot compute statistics.');
end

h_ok = h_final_all(success_flags, :, :);
F_ok = F_final_all(success_flags, :, :);

mc_results = struct();
mc_results.h_mean_final  = squeeze(mean(h_ok, 1));
mc_results.h_std_final   = squeeze(std(h_ok, 0, 1));
mc_results.h_p025_final  = squeeze(prctile(h_ok, 2.5, 1));
mc_results.h_p975_final  = squeeze(prctile(h_ok, 97.5, 1));

mc_results.F_mean_final  = squeeze(mean(F_ok, 1));
mc_results.F_std_final   = squeeze(std(F_ok, 0, 1));

% Outlet Q statistics
if ~isempty(Q_outlet_all)
    Q_ok = Q_outlet_all(success_flags, :);
    mc_results.t_common       = t_common;
    mc_results.Q_outlet_mean  = mean(Q_ok, 1, 'omitnan')';
    mc_results.Q_outlet_std   = std(Q_ok, 0, 1, 'omitnan')';
    mc_results.Q_outlet_p025  = prctile(Q_ok, 2.5, 1)';
    mc_results.Q_outlet_p975  = prctile(Q_ok, 97.5, 1)';

    if save_traj
        mc_results.Q_outlet_ensemble = Q_ok;
    end
else
    mc_results.t_common = [];
    mc_results.Q_outlet_mean = [];
    mc_results.Q_outlet_std = [];
end

mc_results.param_samples  = param_samples;
mc_results.param_names    = param_names;
mc_results.n_samples      = N;
mc_results.n_successful   = n_ok;
mc_results.solver_name    = solver_name;
mc_results.random_seed    = random_seed;
mc_results.sampling_method = sampling_method;
mc_results.computation_time = toc(tic_total);
mc_results.convergence    = convergence;

fprintf('  MC validation complete: %.1f s, %d/%d successful\n', ...
    mc_results.computation_time, n_ok, N);
fprintf('  Mean σ_h(final) = %.4f mm\n', mean(mc_results.h_std_final(:), 'omitnan') * 1000);

end % run_monte_carlo_validation


%% ========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function [names, base_vals, sigma_vals] = identify_active_params(sim_params, pse_config)
% Identify which parameters have non-zero CV in this PSE config

names = {};
base_vals = [];
sigma_vals = [];

% Rainfall CV
if isfield(pse_config, 'rainfall_cv') && pse_config.rainfall_cv > 0
    % Use a representative intensity [m/s]
    R_rep = 50 / 1000 / 3600;  % 50 mm/hr default
    names{end+1} = 'rainfall';
    base_vals(end+1) = R_rep;
    sigma_vals(end+1) = R_rep * pse_config.rainfall_cv;
end

% Manning's n CV
if isfield(pse_config, 'manning_n_cv') && pse_config.manning_n_cv > 0
    if isfield(sim_params, 'n_manning')
        n_vals = sim_params.n_manning(~isnan(sim_params.n_manning));
        n_mean = mean(n_vals(:));
    else
        n_mean = 0.03;
    end
    names{end+1} = 'manning_n';
    base_vals(end+1) = n_mean;
    sigma_vals(end+1) = n_mean * pse_config.manning_n_cv;
end

% Ks CV
if isfield(pse_config, 'Ks_cv') && pse_config.Ks_cv > 0 && isfield(sim_params, 'Ks')
    Ks_vals = sim_params.Ks(~isnan(sim_params.Ks));
    Ks_mean = mean(Ks_vals(:));
    names{end+1} = 'Ks';
    base_vals(end+1) = Ks_mean;
    sigma_vals(end+1) = Ks_mean * pse_config.Ks_cv;
end

% Psi CV (field name: psi_f in domain setup)
if isfield(pse_config, 'psi_cv') && pse_config.psi_cv > 0 && isfield(sim_params, 'psi_f')
    if isnumeric(sim_params.psi_f) && numel(sim_params.psi_f) > 1
        psi_vals = sim_params.psi_f(~isnan(sim_params.psi_f));
        psi_mean = mean(psi_vals(:));
    else
        psi_mean = sim_params.psi_f;
    end
    names{end+1} = 'psi';
    base_vals(end+1) = psi_mean;
    sigma_vals(end+1) = psi_mean * pse_config.psi_cv;
end

% Eta (effective porosity, field name: eta in domain setup)
if isfield(pse_config, 'eta_cv') && pse_config.eta_cv > 0 && isfield(sim_params, 'eta')
    if isnumeric(sim_params.eta) && numel(sim_params.eta) > 1
        eta_vals = sim_params.eta(~isnan(sim_params.eta));
        eta_mean = mean(eta_vals(:));
    else
        eta_mean = sim_params.eta;
    end
    names{end+1} = 'eta';
    base_vals(end+1) = eta_mean;
    sigma_vals(end+1) = eta_mean * pse_config.eta_cv;
end

end


function sp = perturb_params(sim_params, param_names, sample_vals)
% PERTURB_PARAMS  Apply sampled parameter values to a copy of sim_params

sp = sim_params;

for pp = 1:numel(param_names)
    name = param_names{pp};
    val = sample_vals(pp);

    switch name
        case 'manning_n'
            % Scale Manning's n uniformly (multiply all cells by ratio)
            if isfield(sp, 'n_manning')
                n_base = sp.n_manning;
                n_mean = mean(n_base(~isnan(n_base)));
                if n_mean > 0
                    sp.n_manning = n_base * (val / n_mean);
                end
            end

        case 'Ks'
            % Scale Ks field
            if isfield(sp, 'Ks')
                Ks_base = sp.Ks;
                Ks_mean = mean(Ks_base(~isnan(Ks_base)));
                if Ks_mean > 0
                    sp.Ks = Ks_base * (val / Ks_mean);
                end
            end

        case 'psi'
            % Field name: psi_f in domain setup
            if isfield(sp, 'psi_f')
                if numel(sp.psi_f) > 1
                    psi_base = sp.psi_f;
                    psi_mean = mean(psi_base(~isnan(psi_base)));
                    if psi_mean > 0
                        sp.psi_f = psi_base * (val / psi_mean);
                    end
                else
                    sp.psi_f = val;
                end
            end

        case 'eta'
            % Field name: eta in domain setup (not delta_theta)
            if isfield(sp, 'eta')
                if numel(sp.eta) > 1
                    eta_base = sp.eta;
                    eta_mean = mean(eta_base(~isnan(eta_base)));
                    if eta_mean > 0
                        sp.eta = eta_base * (val / eta_mean);
                    end
                else
                    sp.eta = val;
                end
            end

        case 'rainfall'
            % Rainfall handled separately in perturb_rainfall
    end
end
end


function R_out = perturb_rainfall(R_intensity, pse_config, param_names, sample_vals)
% PERTURB_RAINFALL  Scale rainfall by sampled multiplier

R_out = R_intensity;

rain_idx = find(strcmp(param_names, 'rainfall'), 1);
if isempty(rain_idx)
    return;  % Rainfall not being perturbed
end

% Compute scale factor from base intensity
R_rep = 50 / 1000 / 3600;  % Same reference as identify_active_params
rain_val = sample_vals(rain_idx);
scale = rain_val / R_rep;

if isnumeric(R_out)
    R_out = R_out * scale;
elseif isstruct(R_out)
    if isfield(R_out, 'base')
        R_out.base = R_out.base * scale;
    end
    if isfield(R_out, 'schedule')
        R_out.schedule(:, 3) = R_out.schedule(:, 3) * scale;
    end
    if isfield(R_out, 'rain_matrix')
        R_out.rain_matrix = R_out.rain_matrix * scale;
    end
end
end


function param_dist = get_param_distributions(pse_config, param_names)
% GET_PARAM_DISTRIBUTIONS  Extract distribution types for parameters
% Returns cell array matching param_names with distribution type for each

param_dist = cell(numel(param_names), 1);

% Default distributions (match PSE defaults)
default_dist = struct();
default_dist.rainfall = 'lognormal';
default_dist.manning_n = 'lognormal';
default_dist.Ks = 'lognormal';
default_dist.psi = 'lognormal';
default_dist.eta = 'truncated_gaussian';

% Override with user config if present
if isfield(pse_config, 'param_distributions') && isstruct(pse_config.param_distributions)
    user_dist = pse_config.param_distributions;
    fields = fieldnames(user_dist);
    for ff = 1:numel(fields)
        default_dist.(fields{ff}) = user_dist.(fields{ff});
    end
end

% Map parameter names to distribution types
for pp = 1:numel(param_names)
    name = param_names{pp};
    if isfield(default_dist, name)
        param_dist{pp} = default_dist.(name);
    else
        param_dist{pp} = 'lognormal';  % Fallback
    end
end
end
