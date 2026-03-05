function pse_results = run_pse_with_measurements(sim_params, sim_results, t_out, x_out, pse_config, sensor_config)
% RUN_PSE_WITH_MEASUREMENTS  PSE with measurement matrix
%
% MATHEMATICAL FRAMEWORK:
% =======================
% Extends run_pse_dae to include measurement equations following
%  "Probabilistic State Estimation in Water Networks".
%
% The standard PSE operates on the full state vector x and propagates
% uncertainty from parameters to states via:
%
%   K_xx = A^{-1} K_bb A^{-T}
%
% When sensors are available, we augment the system with measurements:
%
%   y(k) = C * x(k) + v(k)     
%
% where C is a binary selection matrix and v is measurement noise.
% The augmented system becomes:
%
%   A^s = [E; C]                 (Eq. 10: DAE + measurements)
%   b^s = [z; y - v]            (RHS: parameters + observations)
%   K_bb^s = blkdiag(K_zz, K_vv)  (Eq. 11: combined covariance)
%
% PSE with measurements (GLS/BLUE estimator, Corollary 1):
%
%   K_xx = (A^T W A)^{-1}   where W = (K_bb^s)^{-1}
%
% Note: When W = inv(K_bb^s), the sandwich formula from Theorem 1
% collapses to this simpler form (Gauss-Markov optimal estimator).
%
% KEY INSIGHT: Measurements reduce state uncertainty. The more sensors
% (rows in C), the more constrained the states become, and the smaller
% K_xx. This is the probabilistic analog of observability.
%
% OUTPUT UNCERTAINTY:
%   K_yy = C * K_xx * C' + K_vv   (measurement-space covariance)
%
% COMPARISON WITH run_pse_dae:
%   - run_pse_dae: K_xx from A and K_bb only (no C, no measurements)
%   - This function: K_xx from [A; C] and blkdiag(K_bb, K_vv)
%   - When sensor_config has no sensors, this reduces exactly to run_pse_dae
%
% STATE VECTOR: x = [h, Q_l, Q_r, Q_u, Q_d, F] (6 states per cell)
%
% INPUTS:
%   sim_params    - Simulation parameters (from setup)
%   sim_results   - DAE simulation results
%   t_out         - Time vector from simulation
%   x_out         - State trajectory matrix (nt × n_states_total)
%   pse_config    - PSE configuration (CVs, confidence level, etc.)
%   sensor_config - Sensor placement from configure_sensors()
%
% OUTPUTS:
%   pse_results   - Structure with fields:
%     .h_std, .F_std, .Q_*_std  - Standard deviations (same as run_pse_dae)
%     .h_lower, .h_upper, ...   - Confidence intervals
%     .y_std                    - Measurement output std (per sensor, per time)
%     .y_lower, .y_upper        - Measurement confidence intervals
%     .sensor_config            - Copy of sensor configuration
%     .variance_reduction       - Ratio of with/without measurement variance
%
% REFERENCE:
%   Sections III-B and III-C
%
% Author: M. H. Kazma
% Date: February 2026

%% Initialize
fprintf('  Initializing PSE WITH MEASUREMENTS...\n');
fprintf('  System: A^s = [E; C], K_bb^s = blkdiag(K_zz, K_vv)\n');
fprintf('  Formula: K_xx = (A^T A)^{-1} A^T K_bb^s A (A^T A)^{-1}\n');

if nargin < 5 || isempty(pse_config)
    pse_config = struct();
end
if nargin < 6 || isempty(sensor_config)
    error('sensor_config is required. Use configure_sensors() to create it.');
end

% Default PSE config values
if ~isfield(pse_config, 'rainfall_cv'), pse_config.rainfall_cv = 0.20; end
if ~isfield(pse_config, 'rainfall_sigma_abs'), pse_config.rainfall_sigma_abs = 0; end
if ~isfield(pse_config, 'manning_n_cv'), pse_config.manning_n_cv = 0.20; end
if ~isfield(pse_config, 'Ks_cv'), pse_config.Ks_cv = 0.20; end
if ~isfield(pse_config, 'psi_cv'), pse_config.psi_cv = 0.20; end
if ~isfield(pse_config, 'eta_cv'), pse_config.eta_cv = 0.20; end
debug_mode = isfield(pse_config, 'debug') && pse_config.debug;

% Parameter distribution configuration (match run_pse_dae)
param_dist = struct();
param_dist.rainfall = 'lognormal';
param_dist.manning_n = 'lognormal';
param_dist.Ks = 'lognormal';
param_dist.psi = 'lognormal';
param_dist.eta = 'truncated_gaussian';
if isfield(pse_config, 'param_distributions') && isstruct(pse_config.param_distributions)
    param_dist = merge_struct_defaults(param_dist, pse_config.param_distributions);
end

param_bounds = struct();
param_bounds.rainfall = [0, inf];
param_bounds.manning_n = [0, inf];
param_bounds.Ks = [0, inf];
param_bounds.psi = [0, inf];
param_bounds.eta = [0, 1];
if isfield(pse_config, 'param_bounds') && isstruct(pse_config.param_bounds)
    param_bounds = merge_struct_defaults(param_bounds, pse_config.param_bounds);
end

%% Temporal measurement schedule
% Controls WHEN sensors provide data. At measurement times, the augmented
% system A^s = [E; C] is used (over-determined → variance reduction).
% At non-measurement times, only physics A^s = E propagates uncertainty
% (square system → variance grows). This creates the "sawtooth" pattern.
if isfield(pse_config, 'measurement_schedule') && isstruct(pse_config.measurement_schedule)
    meas_schedule = pse_config.measurement_schedule;
    use_temporal_schedule = true;
    fprintf('  Temporal scheduling: %s\n', meas_schedule.description);
else
    meas_schedule = [];
    use_temporal_schedule = false;
end

%% Validate inputs
if ~isfield(sim_params, 'routing_method') || ~strcmp(sim_params.routing_method, 'manning_ca')
    error('run_pse_with_measurements supports only manning_ca routing.');
end
if ~isfield(sim_params, 'n_states_per_cell') || sim_params.n_states_per_cell < 6
    error('run_pse_with_measurements expects 6-state formulation.');
end

%% Grid info
ny = sim_params.grid_size(1);
nx = sim_params.grid_size(2);

use_compressed = isfield(sim_params, 'use_compressed_indexing') && sim_params.use_compressed_indexing;
if use_compressed
    compressed_to_grid = sim_params.compressed_to_grid;
    grid_to_compressed = sim_params.grid_to_compressed;
    n_active = sim_params.valid_cells;
    active_mask = false(ny, nx);
    active_mask(sub2ind([ny, nx], compressed_to_grid(:,1), compressed_to_grid(:,2))) = true;
else
    n_active = nx * ny;
    active_mask = true(ny, nx);
    grid_to_compressed = [];
    compressed_to_grid = [];
end

%% Parameters
Z = sim_params.Z;
dx = sim_params.dx;
dy = sim_params.dy;
cell_area = dx * dy;  % Used for neighbor inflow variance (Fix 3)

if isfield(sim_params, 'n_manning')
    n_manning = sim_params.n_manning;
    if isscalar(n_manning)
        n_manning = n_manning * ones(ny, nx);
    end
else
    error('Manning n not found in sim_params');
end

has_infiltration = isfield(sim_params, 'Ks');
if has_infiltration
    Ks = sim_params.Ks;
    eta = sim_params.eta;
    psi_f = sim_params.psi_f;
    F_min = sim_params.F_min;
end

% Outlet/boundary parameters
outlet_type = 1;
outlet_indices = [];
slope_outlet = 0.01;
outlet_direction = 'right';
if isfield(sim_params, 'outlet_type'), outlet_type = sim_params.outlet_type; end
if isfield(sim_params, 'outlet_indices'), outlet_indices = sim_params.outlet_indices; end
if isfield(sim_params, 'slope_outlet'), slope_outlet = sim_params.slope_outlet; end
if isfield(sim_params, 'outlet_direction'), outlet_direction = sim_params.outlet_direction; end

% Regularization parameters
if isfield(sim_params, 'h_smooth'), h_smooth = sim_params.h_smooth;
else
    if isfield(sim_params, 'slope_max') && sim_params.slope_max > 0.2
        h_smooth = 2e-3;
    else
        h_smooth = 1e-4;
    end
end
if isfield(sim_params, 'h_0'), h_0 = sim_params.h_0; else, h_0 = 0; end
if isfield(sim_params, 'h_flow_threshold'), h_flow_threshold = sim_params.h_flow_threshold; else, h_flow_threshold = 1e-3; end
if isfield(sim_params, 'v_limit'), v_limit = sim_params.v_limit; else, v_limit = 8.0; end

%% Time discretization
nt = length(t_out);

% Optional: resample trajectory to uniform time step
if isfield(pse_config, 'resample_dt_s') && pse_config.resample_dt_s > 0
    dt_resample = pse_config.resample_dt_s;
    t_uniform = (t_out(1):dt_resample:t_out(end))';
    if t_uniform(end) < t_out(end)
        t_uniform = [t_uniform; t_out(end)];
    end
    x_out = interp1(t_out(:), x_out, t_uniform, 'linear', 'extrap');
    t_out = t_uniform;
    nt = length(t_out);
end

dt_vec = diff(t_out);

%% Rainfall extraction (reuse from run_pse_dae, with rain_matrix support)
if ~isfield(sim_results, 'R')
    error('Rainfall data not found in sim_results');
end
R_in = sim_results.R;

use_rain_matrix_pse = false;
R_rain_matrix_grid = [];  % [ny x nx x nt] if rain_matrix mode
R_time = [];
R_pattern_scalar = 1;
R_pattern_grid = [];

if isstruct(R_in) && isfield(R_in, 'rain_matrix')
    % Spatially distributed, time-varying rain matrix [n_cells x n_rain_steps]
    rm = R_in.rain_matrix;
    dt_rain = R_in.dt_rain;
    use_rain_matrix_pse = true;
    R_rain_matrix_grid = zeros(ny, nx, nt);
    for kk = 1:nt
        t_idx = min(floor(t_out(kk) / dt_rain) + 1, size(rm, 2));
        for cc = 1:n_active
            row_c = compressed_to_grid(cc, 1);
            col_c = compressed_to_grid(cc, 2);
            R_rain_matrix_grid(row_c, col_c, kk) = rm(cc, t_idx);
        end
    end
    R_time = squeeze(mean(mean(R_rain_matrix_grid, 1), 2))';
    fprintf('  PSE: rain_matrix mode (%d x %d x %d)\n', ny, nx, nt);
else
    [R_time, R_pattern_scalar, R_pattern_grid] = extract_rainfall(R_in, t_out, nt, ...
        use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
end

if isempty(R_time)
    R_time = ones(1, nt);
end

%% Extract deterministic trajectories
fprintf('  Extracting deterministic trajectories...\n');
n_states_per_cell = sim_params.n_states_per_cell;
[h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj] = ...
    extract_state_trajectories(x_out, t_out, nt, ny, nx, n_states_per_cell, ...
    use_compressed, n_active, compressed_to_grid);

%% Preallocate state variances
var_h = zeros(ny, nx, nt);
var_F = zeros(ny, nx, nt);
var_Ql = zeros(ny, nx, nt);
var_Qr = zeros(ny, nx, nt);
var_Qu = zeros(ny, nx, nt);
var_Qd = zeros(ny, nx, nt);
cov_hF = zeros(ny, nx, nt);  % Fix 1: Track h-F cross-covariance for temporal propagation

% Preallocate measurement output variances
n_sensors_total = sensor_config.n_total;
var_y = zeros(n_sensors_total, nt);  % Per-sensor variance over time

% Build a lookup: for each sensor, which cell and which state index
sensor_cell_map = build_sensor_cell_map(sensor_config, ny, nx);

%% Pre-compute variance scaling factors (match run_pse_dae optimizations)
cv_n = pse_config.manning_n_cv;

% Pre-compute sigma_ln_Ks ONCE (avoid isfield/sqrt/log in inner loop)
if has_infiltration
    if isfield(pse_config, 'sigma_ln_Ks')
        sig_lnKs = pse_config.sigma_ln_Ks;
    else
        sig_lnKs = sqrt(log(1 + pse_config.Ks_cv^2));
    end
else
    sig_lnKs = 0;
end

% PSE computation stride: compute every N timesteps to speed up
if isfield(pse_config, 'stride') && pse_config.stride > 1
    pse_stride = pse_config.stride;
else
    pse_stride = 1;
end

% Use analytical Jacobian unless user forces numerical
use_analytical_jacobian = true;
if isfield(pse_config, 'use_numerical_jacobian') && pse_config.use_numerical_jacobian
    use_analytical_jacobian = false;
end

%% ============================================================
%  MAIN LOOP: PSE with measurement augmentation
%  ============================================================
if use_analytical_jacobian
    fprintf('  Computing K_xx with ANALYTICAL Jacobian + measurements...\n');
else
    fprintf('  Computing K_xx with NUMERICAL Jacobian + measurements...\n');
end
fprintf('  Sensors: %d depth, %d flow, %d infil (total %d)\n', ...
    sensor_config.n_depth, sensor_config.n_flow, sensor_config.n_infil, n_sensors_total);
if use_temporal_schedule
    fprintf('  Temporal schedule: %d/%d timesteps with measurements (%.1f%%)\n', ...
        meas_schedule.n_measurements, meas_schedule.n_total, 100*meas_schedule.fraction);
else
    fprintf('  Temporal schedule: continuous (measurements at every timestep)\n');
end
fprintf('  Stride=%d, %d active cells, %d timesteps\n', pse_stride, n_active, nt);

tic;
wall_start = tic;

active_indices = find(active_mask);
[active_rows, active_cols] = ind2sub([ny, nx], active_indices);

% Track statistics
neg_diag_count = 0;
total_cells = 0;
measured_cells = 0;

% Weight matrix from sensor config
W = sensor_config.weight_matrix;

if debug_mode
    min_rcond = inf;
    variance_reduction_sum = 0;
    variance_reduction_count = 0;
end

% Determine which timesteps to compute PSE for
pse_timesteps = 2:pse_stride:nt;
% Always include the last timestep
if pse_timesteps(end) ~= nt
    pse_timesteps = [pse_timesteps, nt];
end
n_pse_steps = length(pse_timesteps);
pse_step_count = 0;

for ki = 1:n_pse_steps
    k = pse_timesteps(ki);
    % Previous PSE timestep index (for variance retrieval)
    if ki > 1
        k_prev_pse = pse_timesteps(ki-1);
    else
        k_prev_pse = 1;  % initial conditions
    end
    % For stride > 1, accumulate dt across skipped steps
    if pse_stride > 1 && ki > 1
        dt = sum(dt_vec(k_prev_pse:k-1));
    else
        dt = dt_vec(k-1);
    end

    pse_step_count = pse_step_count + 1;
    if mod(pse_step_count, 10) == 0 || pse_step_count == 1
        elapsed = toc(wall_start);
        frac = pse_step_count / n_pse_steps;
        if frac > 0.01
            eta_sec = elapsed / frac * (1 - frac);
            fprintf('    PSE step %d/%d (t=%d/%d) [%.0fs elapsed, ~%.0fs remaining]\n', ...
                pse_step_count, n_pse_steps, k, nt, elapsed, eta_sec);
        else
            fprintf('    PSE step %d/%d (t=%d/%d)\n', pse_step_count, n_pse_steps, k, nt);
        end
    end

    for a = 1:numel(active_rows)
        row = active_rows(a);
        col = active_cols(a);

        % Quick dry cell check FIRST — avoid expensive operations for dry cells
        h_curr = h_traj(row, col, k);
        dry_threshold = 1e-8;
        dry_decay_rain = 0.5;
        dry_decay_norain = 0.1;
        if isfield(pse_config, 'dry_cell_threshold'), dry_threshold = pse_config.dry_cell_threshold; end
        if isfield(pse_config, 'dry_cell_decay_rain'), dry_decay_rain = pse_config.dry_cell_decay_rain; end
        if isfield(pse_config, 'dry_cell_decay_norain'), dry_decay_norain = pse_config.dry_cell_decay_norain; end

        if h_curr < dry_threshold
            % --- DRY CELL FAST PATH (matches run_pse_dae) ---
            var_h_prev = var_h(row, col, k_prev_pse);
            var_F_prev = var_F(row, col, k_prev_pse);

            % Get rainfall for this cell
            if use_rain_matrix_pse
                R_cell = R_rain_matrix_grid(row, col, k);
            elseif ~isempty(R_pattern_grid)
                R_cell = R_pattern_grid(row, col) * R_time(k);
            else
                R_cell = R_pattern_scalar * R_time(k);
            end

            var_R_dry = (dt * pse_config.rainfall_cv * R_cell)^2 + ...
                        (dt * pse_config.rainfall_sigma_abs)^2;

            % F variance: pre-ponding, F integrates rainfall
            if has_infiltration
                Ks_cell = get_cell_param(Ks, row, col);
                var_Ks_contrib = (Ks_cell * sig_lnKs * dt)^2;
                var_F(row, col, k) = var_F_prev + var_R_dry + var_Ks_contrib;

                if R_cell > 0
                    var_excess = var_R_dry + var_Ks_contrib;
                    var_h(row, col, k) = max(var_h_prev * dry_decay_rain, var_excess);
                else
                    var_h(row, col, k) = var_h_prev * dry_decay_norain;
                end
            else
                var_F(row, col, k) = var_F_prev + var_R_dry;
                if R_cell > 0
                    var_h(row, col, k) = var_R_dry;
                else
                    var_h(row, col, k) = var_h_prev * dry_decay_norain;
                end
            end
            % Q: zero for dry cells (already initialized to 0)
            cov_hF(row, col, k) = 0;  % Fix 1: Reset cross-covariance for dry cells
            continue;
        end

        % --- WET CELL PATH ---
        F_curr = F_traj(row, col, k);
        Q_l = Q_l_traj(row, col, k);
        Q_r = Q_r_traj(row, col, k);
        Q_u = Q_u_traj(row, col, k);
        Q_d = Q_d_traj(row, col, k);

        y0 = [h_curr; Q_l; Q_r; Q_u; Q_d; F_curr];

        xdot_h = (h_curr - h_traj(row, col, k_prev_pse)) / dt;
        xdot_F = (F_curr - F_traj(row, col, k_prev_pse)) / dt;

        n_cell = n_manning(row, col);
        if use_rain_matrix_pse
            R_cell = R_rain_matrix_grid(row, col, k);
        elseif ~isempty(R_pattern_grid)
            R_cell = R_pattern_grid(row, col) * R_time(k);
        else
            R_cell = R_pattern_scalar * R_time(k);
        end
        if has_infiltration
            Ks_cell = get_cell_param(Ks, row, col);
            psi_cell = get_cell_param(psi_f, row, col);
            eta_cell = get_cell_param(eta, row, col);
        else
            Ks_cell = 0; psi_cell = 0; eta_cell = 0;
        end

        % Previous variances for temporal propagation
        var_h_prev = var_h(row, col, k_prev_pse);
        var_F_prev = var_F(row, col, k_prev_pse);
        cov_hF_prev = cov_hF(row, col, k_prev_pse);  % Fix 1: Previous h-F cross-covariance

        %% Build 6x6 DAE Jacobian (analytical or numerical)
        if use_analytical_jacobian
            A_dae = build_dae_jacobian_6x6_analytical_local(y0, dt, ...
                n_cell, Ks_cell, psi_cell, eta_cell, ...
                dx, dy, h_0, h_smooth, h_flow_threshold, v_limit, F_min, has_infiltration);
        else
            A_dae = build_dae_jacobian_6x6(y0, row, col, k, dt, xdot_h, xdot_F, ...
                R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
                h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, ...
                Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, ...
                active_mask, F_min, has_infiltration, ...
                outlet_type, outlet_indices, outlet_direction, slope_outlet, ...
                use_compressed, grid_to_compressed);
        end

        %% Build K_bb (6x6 covariance of DAE RHS) — same as run_pse_dae
        K_bb = build_kbb_dae(h_curr, F_curr, Q_l, Q_r, Q_u, Q_d, ...
            var_h_prev, var_F_prev, cov_hF_prev, dt, R_cell, n_cell, ...
            Ks_cell, psi_cell, eta_cell, F_min, has_infiltration, pse_config, sig_lnKs);

        %% Fix 3: Neighbor inflow uncertainty
        % In the paper's global system, neighbor Q uncertainty propagates automatically.
        % For cell-local PSE, we inject variance of inflow Q from neighbors into K_bb(1,1).
        var_Q_in = 0;
        if col > 1 && active_mask(row, col-1)
            var_Q_in = var_Q_in + var_Qr(row, col-1, k_prev_pse);
        end
        if col < nx && active_mask(row, col+1)
            var_Q_in = var_Q_in + var_Ql(row, col+1, k_prev_pse);
        end
        if row > 1 && active_mask(row-1, col)
            var_Q_in = var_Q_in + var_Qd(row-1, col, k_prev_pse);
        end
        if row < ny && active_mask(row+1, col)
            var_Q_in = var_Q_in + var_Qu(row+1, col, k_prev_pse);
        end
        K_bb(1,1) = K_bb(1,1) + dt^2 * var_Q_in / cell_area^2;

        %% Build measurement matrix C and noise K_vv for this cell
        % Temporal scheduling: only include measurements at scheduled times
        if use_temporal_schedule && ~meas_schedule.is_measurement(k)
            % Non-measurement timestep: pure physics propagation
            % A^s = E (square), K_xx = A^{-1} K_bb A^{-T}
            C_local = [];
            K_vv_local = [];
        else
            % Measurement timestep (or no scheduling): augmented system
            % A^s = [E; C] (over-determined), K_xx = (A^T W A)^{-1}
            [C_local, K_vv_local, ~] = build_measurement_matrix(row, col, sensor_config);
        end

        %% Compute K_xx with measurement augmentation
        [K_xx, K_yy_local, rc, had_neg_diag] = compute_kxx_measured(...
            A_dae, K_bb, C_local, K_vv_local, W);

        total_cells = total_cells + 1;
        if had_neg_diag
            neg_diag_count = neg_diag_count + 1;
        end

        %% Store measurement output variances
        if ~isempty(C_local) && ~isempty(K_yy_local)
            measured_cells = measured_cells + 1;
            var_y = store_sensor_variances(var_y, K_yy_local, row, col, k, sensor_config, sensor_cell_map);

            if debug_mode
                K_xx_nomeas = compute_kxx_measured(A_dae, K_bb, [], [], []);
                for si = 1:size(C_local, 1)
                    state_idx = find(C_local(si, :));
                    if ~isempty(state_idx)
                        var_meas = max(0, K_xx(state_idx, state_idx));
                        var_nomeas = max(0, K_xx_nomeas(state_idx, state_idx));
                        if var_nomeas > 1e-20
                            reduction = var_meas / var_nomeas;
                            variance_reduction_sum = variance_reduction_sum + reduction;
                            variance_reduction_count = variance_reduction_count + 1;
                        end
                    end
                end
            end
        end

        if debug_mode && rc < min_rcond
            min_rcond = rc;
        end

        %% Extract state variances from diagonal of K_xx
        var_h(row, col, k) = max(0, K_xx(1,1));
        var_Ql(row, col, k) = max(0, K_xx(2,2));
        var_Qr(row, col, k) = max(0, K_xx(3,3));
        var_Qu(row, col, k) = max(0, K_xx(4,4));
        var_Qd(row, col, k) = max(0, K_xx(5,5));
        var_F(row, col, k) = max(0, K_xx(6,6));
        cov_hF(row, col, k) = K_xx(1,6);  % Fix 1: Store h-F cross-covariance for next timestep
    end
end

% For stride > 1, interpolate skipped timesteps
if pse_stride > 1
    stride_interp = 'linear';
    if isfield(pse_config, 'stride_interp'), stride_interp = pse_config.stride_interp; end
    fprintf('  Interpolating %d skipped timesteps (stride=%d, method=%s)...\n', ...
        nt - n_pse_steps - 1, pse_stride, stride_interp);

    if strcmp(stride_interp, 'hold')
        % Zero-order hold: copy previous computed value to skipped steps
        for ki = 2:n_pse_steps
            k_curr = pse_timesteps(ki);
            k_prev = pse_timesteps(ki-1);
            if k_curr - k_prev > 1
                for kk = (k_prev+1):(k_curr-1)
                    var_h(:,:,kk) = var_h(:,:,k_prev);
                    var_F(:,:,kk) = var_F(:,:,k_prev);
                    var_Ql(:,:,kk) = var_Ql(:,:,k_prev);
                    var_Qr(:,:,kk) = var_Qr(:,:,k_prev);
                    var_Qu(:,:,kk) = var_Qu(:,:,k_prev);
                    var_Qd(:,:,kk) = var_Qd(:,:,k_prev);
                    cov_hF(:,:,kk) = cov_hF(:,:,k_prev);
                end
            end
        end
    elseif strcmp(stride_interp, 'pchip')
        % PCHIP interpolation per-cell (shape-preserving, C^1)
        t_pse = t_out(pse_timesteps);
        for a = 1:numel(active_rows)
            r = active_rows(a); c = active_cols(a);
            vals_h  = squeeze(var_h(r,c,pse_timesteps));
            vals_F  = squeeze(var_F(r,c,pse_timesteps));
            vals_Ql = squeeze(var_Ql(r,c,pse_timesteps));
            vals_Qr = squeeze(var_Qr(r,c,pse_timesteps));
            vals_Qu = squeeze(var_Qu(r,c,pse_timesteps));
            vals_Qd = squeeze(var_Qd(r,c,pse_timesteps));
            vals_cov = squeeze(cov_hF(r,c,pse_timesteps));
            for ki = 2:n_pse_steps
                k_curr = pse_timesteps(ki);
                k_prev = pse_timesteps(ki-1);
                if k_curr - k_prev > 1
                    kk_range = (k_prev+1):(k_curr-1);
                    t_interp = t_out(kk_range);
                    var_h(r,c,kk_range) = max(0, interp1(t_pse, vals_h, t_interp, 'pchip'));
                    var_F(r,c,kk_range) = max(0, interp1(t_pse, vals_F, t_interp, 'pchip'));
                    var_Ql(r,c,kk_range) = max(0, interp1(t_pse, vals_Ql, t_interp, 'pchip'));
                    var_Qr(r,c,kk_range) = max(0, interp1(t_pse, vals_Qr, t_interp, 'pchip'));
                    var_Qu(r,c,kk_range) = max(0, interp1(t_pse, vals_Qu, t_interp, 'pchip'));
                    var_Qd(r,c,kk_range) = max(0, interp1(t_pse, vals_Qd, t_interp, 'pchip'));
                    cov_hF(r,c,kk_range) = interp1(t_pse, vals_cov, t_interp, 'pchip');
                end
            end
        end
    else
        % Linear interpolation (default)
        for ki = 2:n_pse_steps
            k_curr = pse_timesteps(ki);
            k_prev = pse_timesteps(ki-1);
            if k_curr - k_prev > 1
                for kk = (k_prev+1):(k_curr-1)
                    alpha = (kk - k_prev) / (k_curr - k_prev);
                    var_h(:,:,kk) = (1-alpha) * var_h(:,:,k_prev) + alpha * var_h(:,:,k_curr);
                    var_F(:,:,kk) = (1-alpha) * var_F(:,:,k_prev) + alpha * var_F(:,:,k_curr);
                    var_Ql(:,:,kk) = (1-alpha) * var_Ql(:,:,k_prev) + alpha * var_Ql(:,:,k_curr);
                    var_Qr(:,:,kk) = (1-alpha) * var_Qr(:,:,k_prev) + alpha * var_Qr(:,:,k_curr);
                    var_Qu(:,:,kk) = (1-alpha) * var_Qu(:,:,k_prev) + alpha * var_Qu(:,:,k_curr);
                    var_Qd(:,:,kk) = (1-alpha) * var_Qd(:,:,k_prev) + alpha * var_Qd(:,:,k_curr);
                    cov_hF(:,:,kk) = (1-alpha) * cov_hF(:,:,k_prev) + alpha * cov_hF(:,:,k_curr);  % Fix 1
                end
            end
        end
    end
end

pse_wall_time = toc(wall_start);
computation_time = toc;

%% Summary and Diagnostics
fprintf('  Variance propagation complete in %.1f seconds.\n', pse_wall_time);
if use_analytical_jacobian
    fprintf('  - Method: analytical Jacobian + measurements, stride=%d\n', pse_stride);
else
    fprintf('  - Method: numerical Jacobian + measurements, stride=%d\n', pse_stride);
end
fprintf('  - h, F: Temporal propagation + parameter injection\n');
fprintf('  - Q: Parameter injection only (algebraic, no temporal memory)\n');
fprintf('  Measured cells processed: %d\n', measured_cells);

% PSE coverage diagnostics
var_h_final = var_h(:,:,nt);
var_F_final = var_F(:,:,nt);
var_Ql_final = var_Ql(:,:,nt);
n_h_nonzero = nnz(var_h_final(active_mask) > 0);
n_F_nonzero = nnz(var_F_final(active_mask) > 0);
n_Ql_nonzero = nnz(var_Ql_final(active_mask) > 0);
fprintf('  PSE coverage at final timestep:\n');
fprintf('    var_h > 0: %d/%d cells (%.1f%%)\n', n_h_nonzero, n_active, 100*n_h_nonzero/max(n_active,1));
fprintf('    var_F > 0: %d/%d cells (%.1f%%)\n', n_F_nonzero, n_active, 100*n_F_nonzero/max(n_active,1));
fprintf('    var_Ql > 0: %d/%d cells (%.1f%%)\n', n_Ql_nonzero, n_active, 100*n_Ql_nonzero/max(n_active,1));

max_h_std_ever = max(sqrt(var_h(:)));
max_F_std_ever = max(sqrt(var_F(:)));
max_Ql_std_ever = max(sqrt(var_Ql(:)));
fprintf('    max h_std: %.6f m\n', max_h_std_ever);
fprintf('    max F_std: %.6f m\n', max_F_std_ever);
fprintf('    max Q_l_std: %.6f m³/s\n', max_Ql_std_ever);

%% Package results
h_std = sqrt(var_h);
F_std = sqrt(var_F);
Q_l_std = sqrt(var_Ql);
Q_r_std = sqrt(var_Qr);
Q_u_std = sqrt(var_Qu);
Q_d_std = sqrt(var_Qd);

pse_results = struct();

% Standard deviations (from PSE covariance diagonal)
pse_results.h_std = h_std;
pse_results.F_std = F_std;
pse_results.Q_l_std = Q_l_std;
pse_results.Q_r_std = Q_r_std;
pse_results.Q_u_std = Q_u_std;
pse_results.Q_d_std = Q_d_std;

% Net flow: Q_net = Q_l + Q_r + Q_u + Q_d (sum of all directional flows)
% Assuming independent uncertainties: var_net = sum of variances
Q_net_std = sqrt(var_Ql + var_Qr + var_Qu + var_Qd);
pse_results.Q_net_std = Q_net_std;

% Full variance arrays (needed for upstream aggregation in plotting)
pse_results.var_h = var_h;
pse_results.var_F = var_F;
pse_results.var_Ql = var_Ql;
pse_results.var_Qr = var_Qr;
pse_results.var_Qu = var_Qu;
pse_results.var_Qd = var_Qd;

% Measurement output uncertainty
pse_results.var_y = var_y;                           % (n_sensors × nt)
pse_results.y_std = sqrt(max(0, var_y));             % (n_sensors × nt)

% Sensor configuration (for downstream analysis)
pse_results.sensor_config = sensor_config;

%% Confidence intervals
dist_type = 'gaussian';
conf_level = 0.95;
if isfield(pse_config, 'distribution_type'), dist_type = pse_config.distribution_type; end
if isfield(pse_config, 'confidence_level'), conf_level = pse_config.confidence_level; end

fprintf('  Computing %s confidence intervals (%.0f%% level)...\n', dist_type, conf_level*100);

ci_opts = {};
if isfield(pse_config, 'min_depth_for_gaussian')
    ci_opts = [ci_opts, {'MinDepth', pse_config.min_depth_for_gaussian}];
end
if isfield(pse_config, 'truncation_lower')
    ci_opts = [ci_opts, {'TruncLower', pse_config.truncation_lower}];
end

[pse_results.h_lower, pse_results.h_upper] = compute_confidence_intervals(...
    h_traj, h_std, dist_type, conf_level, ci_opts{:});
[pse_results.F_lower, pse_results.F_upper] = compute_confidence_intervals(...
    F_traj, F_std, dist_type, conf_level, ci_opts{:});
[pse_results.Q_l_lower, pse_results.Q_l_upper] = compute_confidence_intervals(...
    Q_l_traj, Q_l_std, 'gaussian', conf_level);
[pse_results.Q_r_lower, pse_results.Q_r_upper] = compute_confidence_intervals(...
    Q_r_traj, Q_r_std, 'gaussian', conf_level);
[pse_results.Q_u_lower, pse_results.Q_u_upper] = compute_confidence_intervals(...
    Q_u_traj, Q_u_std, 'gaussian', conf_level);
[pse_results.Q_d_lower, pse_results.Q_d_upper] = compute_confidence_intervals(...
    Q_d_traj, Q_d_std, 'gaussian', conf_level);

% Net flow CI: Q_net_traj = Q_l + Q_r + Q_u + Q_d
Q_net_traj = Q_l_traj + Q_r_traj + Q_u_traj + Q_d_traj;
[pse_results.Q_net_lower, pse_results.Q_net_upper] = compute_confidence_intervals(...
    Q_net_traj, Q_net_std, 'gaussian', conf_level);

% Measurement output CIs (per sensor)
if n_sensors_total > 0
    % Build measurement mean trajectory from deterministic states
    y_mean = build_measurement_trajectory(h_traj, Q_l_traj, Q_r_traj, ...
        Q_u_traj, Q_d_traj, F_traj, sensor_config, sensor_cell_map, nt);
    pse_results.y_mean = y_mean;
    
    if strcmp(dist_type, 'gaussian')
        z_alpha = norminv(0.5 + conf_level/2);
    else
        z_alpha = 1.96;
    end
    pse_results.y_lower = y_mean - z_alpha * pse_results.y_std;
    pse_results.y_upper = y_mean + z_alpha * pse_results.y_std;
end

%% Variance reduction analysis
if debug_mode && variance_reduction_count > 0
    avg_reduction = variance_reduction_sum / variance_reduction_count;
    fprintf('  Average variance reduction at measured cells: %.1f%%\n', ...
        (1 - avg_reduction) * 100);
    pse_results.avg_variance_reduction = avg_reduction;
end

%% Metadata
pse_results.distribution_type = dist_type;
pse_results.confidence_level = conf_level;
pse_results.param_distributions = param_dist;
pse_results.param_bounds = param_bounds;
pse_results.t = t_out;
pse_results.computation_time = computation_time;
pse_results.method = sprintf('PSE with Measurements (A^s=[E;C]) - %s distribution', dist_type);
pse_results.n_cells = n_active;
pse_results.n_timesteps = nt;
pse_results.grid_size = [ny, nx];
pse_results.n_sensors = n_sensors_total;

fprintf('  Eigenvalue cleanup: %d/%d cells (%.1f%% had negative diagonals)\n', ...
    neg_diag_count, total_cells, 100*neg_diag_count/max(total_cells,1));
fprintf('  PSE with measurements completed in %.2f seconds\n', computation_time);
if debug_mode
    fprintf('  Debug: min rcond = %.2e\n', min_rcond);
end

end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function val = get_cell_param(param, row, col)
if numel(param) > 1
    val = param(row, col);
else
    val = param;
end
end

function K_bb = build_kbb_dae(h_curr, F_curr, Q_l, Q_r, Q_u, Q_d, ...
    var_h_prev, var_F_prev, cov_hF_prev, dt, R_cell, n_cell, ...
    Ks_cell, psi_cell, eta_cell, F_min, has_infiltration, pse_config, sig_lnKs)
% Build the 6x6 covariance of the DAE RHS vector b
% Same logic as in run_pse_dae.m
% sig_lnKs: pre-computed sigma of ln(Ks) for performance
% cov_hF_prev: h-F cross-covariance from previous timestep (Fix 1)

K_bb = zeros(6, 6);

% Temporal propagation: C * K_xx^n * C' where C = diag([1,0,0,0,0,1])
K_bb(1,1) = var_h_prev;
K_bb(6,6) = var_F_prev;
% Fix 1: Include h-F cross-covariance from previous timestep
K_bb(1,6) = cov_hF_prev;
K_bb(6,1) = cov_hF_prev;

% Rainfall uncertainty
% ALWAYS include absolute floor (sigma_abs) even when R=0:
var_R_cell = (dt * pse_config.rainfall_cv * R_cell)^2 + (dt * pse_config.rainfall_sigma_abs)^2;
K_bb(1,1) = K_bb(1,1) + var_R_cell;

% Infiltration uncertainty
if has_infiltration && F_curr > F_min
    f_infil = Ks_cell * (1 + eta_cell * (h_curr + psi_cell) / max(F_curr, F_min));
    
    % Use pre-computed sig_lnKs (log-transform for Ks)
    sigma_ln_Ks = sig_lnKs;
    df_dlnKs = f_infil;
    var_f_from_Ks = df_dlnKs^2 * sigma_ln_Ks^2;
    
    df_dpsi = Ks_cell * eta_cell / max(F_curr, F_min);
    df_deta = Ks_cell * (h_curr + psi_cell) / max(F_curr, F_min);
    var_psi = (pse_config.psi_cv * psi_cell)^2;
    var_eta = (pse_config.eta_cv * eta_cell)^2;
    
    var_f_rate = var_f_from_Ks + df_dpsi^2 * var_psi + df_deta^2 * var_eta;
    var_f_infil = dt^2 * var_f_rate;
else
    var_f_infil = 0;
end

K_bb(1,1) = K_bb(1,1) + var_f_infil;
K_bb(6,6) = K_bb(6,6) + var_f_infil;
K_bb(1,6) = K_bb(1,6) - var_f_infil;   % Add infiltration anti-correlation to temporal cross-cov
K_bb(6,1) = K_bb(1,6);                 % Symmetric

% Manning n affects all Q equations
cv_n = pse_config.manning_n_cv;
var_n = (cv_n * n_cell)^2;
dQl_dn = -Q_l / n_cell;
dQr_dn = -Q_r / n_cell;
dQu_dn = -Q_u / n_cell;
dQd_dn = -Q_d / n_cell;

K_bb(2,2) = dQl_dn^2 * var_n;
K_bb(3,3) = dQr_dn^2 * var_n;
K_bb(4,4) = dQu_dn^2 * var_n;
K_bb(5,5) = dQd_dn^2 * var_n;

K_bb(2,3) = dQl_dn * dQr_dn * var_n;  K_bb(3,2) = K_bb(2,3);
K_bb(2,4) = dQl_dn * dQu_dn * var_n;  K_bb(4,2) = K_bb(2,4);
K_bb(2,5) = dQl_dn * dQd_dn * var_n;  K_bb(5,2) = K_bb(2,5);
K_bb(3,4) = dQr_dn * dQu_dn * var_n;  K_bb(4,3) = K_bb(3,4);
K_bb(3,5) = dQr_dn * dQd_dn * var_n;  K_bb(5,3) = K_bb(3,5);
K_bb(4,5) = dQu_dn * dQd_dn * var_n;  K_bb(5,4) = K_bb(4,5);
end

function [R_time, R_pattern_scalar, R_pattern_grid] = extract_rainfall(R_in, t_out, nt, ...
    use_compressed, n_active, compressed_to_grid, ny, nx, active_mask)
% Extract rainfall time series and spatial pattern (mirrors run_pse_dae logic)
R_time = [];
R_pattern_scalar = 1;
R_pattern_grid = [];

if isstruct(R_in)
    if isfield(R_in, 'schedule')
        schedule = R_in.schedule;
        R_time = zeros(1, nt);
        for kk = 1:nt
            R_time(kk) = get_rainfall_from_schedule(schedule, t_out(kk));
        end
        if isfield(R_in, 'pattern') && ~isempty(R_in.pattern)
            [R_pattern_scalar, R_pattern_grid] = build_rain_pattern(R_in.pattern, ...
                use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
        end
    elseif isfield(R_in, 'base') && isfield(R_in, 'duration')
        base = R_in.base;
        [R_pattern_scalar, R_pattern_grid] = build_rain_pattern(base, ...
            use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
        R_time = ones(1, nt);
        R_time(t_out > R_in.duration) = 0;
    elseif isfield(R_in, 'series') && isfield(R_in, 't')
        R_time = interp1(R_in.t(:)', R_in.series(:)', t_out(:)', 'previous', 'extrap');
    else
        R_pattern_scalar = 1;
        R_time = ones(1, nt);
    end
elseif isscalar(R_in)
    R_pattern_scalar = R_in;
    R_time = ones(1, nt);
elseif isvector(R_in) && length(R_in) == nt
    R_time = R_in(:)';
    R_pattern_scalar = 1;
else
    [R_pattern_scalar, R_pattern_grid] = build_rain_pattern(R_in, ...
        use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
    R_time = ones(1, nt);
end

if isempty(R_time)
    R_time = ones(1, nt);
end
end

function [pattern_scalar, pattern_grid] = build_rain_pattern(base, use_compressed, n_active, compressed_to_grid, ny, nx, active_mask)
pattern_grid = [];
if isscalar(base)
    pattern_scalar = base;
    return;
end
pattern_scalar = 1;
if isvector(base)
    if use_compressed && length(base) == n_active
        pattern_grid = zeros(ny, nx);
        for c = 1:n_active
            pattern_grid(compressed_to_grid(c,1), compressed_to_grid(c,2)) = base(c);
        end
    elseif length(base) == nx*ny
        pattern_grid = reshape(base, [ny, nx]);
    else
        pattern_grid = zeros(ny, nx);
        pattern_grid(active_mask) = base(1);
    end
elseif size(base,1) == ny && size(base,2) == nx
    pattern_grid = base;
else
    pattern_grid = zeros(ny, nx);
    pattern_grid(active_mask) = base(1);
end
end

function [h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj] = ...
    extract_state_trajectories(x_out, t_out, nt, ny, nx, n_states_per_cell, ...
    use_compressed, n_active, compressed_to_grid)
% Extract per-state trajectory arrays from flat x_out

h_traj = zeros(ny, nx, nt);
F_traj = zeros(ny, nx, nt);
Q_l_traj = zeros(ny, nx, nt);
Q_r_traj = zeros(ny, nx, nt);
Q_u_traj = zeros(ny, nx, nt);
Q_d_traj = zeros(ny, nx, nt);

for k = 1:nt
    x_k = x_out(k, :)';
    if use_compressed
        for c = 1:n_active
            i = compressed_to_grid(c,1);
            j = compressed_to_grid(c,2);
            idx_base = (c - 1) * n_states_per_cell;
            h_traj(i,j,k) = x_k(idx_base + 1);
            Q_l_traj(i,j,k) = x_k(idx_base + 2);
            Q_r_traj(i,j,k) = x_k(idx_base + 3);
            Q_u_traj(i,j,k) = x_k(idx_base + 4);
            Q_d_traj(i,j,k) = x_k(idx_base + 5);
            F_traj(i,j,k) = x_k(idx_base + 6);
        end
    else
        for idx = 1:(nx*ny)
            i = ceil(idx / nx);
            j = mod(idx - 1, nx) + 1;
            idx_base = (idx - 1) * n_states_per_cell;
            h_traj(i,j,k) = x_k(idx_base + 1);
            Q_l_traj(i,j,k) = x_k(idx_base + 2);
            Q_r_traj(i,j,k) = x_k(idx_base + 3);
            Q_u_traj(i,j,k) = x_k(idx_base + 4);
            Q_d_traj(i,j,k) = x_k(idx_base + 5);
            F_traj(i,j,k) = x_k(idx_base + 6);
        end
    end
end
end

function sensor_cell_map = build_sensor_cell_map(sensor_config, ny, nx)
% Build a lookup structure mapping each global sensor index to its cell
% and state index for efficient variance storage

n_total = sensor_config.n_total;
sensor_cell_map = struct();
sensor_cell_map.row = zeros(n_total, 1);
sensor_cell_map.col = zeros(n_total, 1);
sensor_cell_map.state_idx = zeros(n_total, 1);  % 1=h, 2=Ql, ..., 6=F
sensor_cell_map.type = cell(n_total, 1);

dir_to_state = [2, 3, 4, 5];  % left→Ql, right→Qr, up→Qu, down→Qd

idx = 0;

% Depth sensors
for si = 1:sensor_config.n_depth
    idx = idx + 1;
    sensor_cell_map.row(idx) = sensor_config.depth_sensors(si, 1);
    sensor_cell_map.col(idx) = sensor_config.depth_sensors(si, 2);
    sensor_cell_map.state_idx(idx) = 1;
    sensor_cell_map.type{idx} = 'depth';
end

% Flow sensors
for si = 1:sensor_config.n_flow
    idx = idx + 1;
    sensor_cell_map.row(idx) = sensor_config.flow_sensors(si, 1);
    sensor_cell_map.col(idx) = sensor_config.flow_sensors(si, 2);
    dir = sensor_config.flow_sensors(si, 3);
    sensor_cell_map.state_idx(idx) = dir_to_state(dir);
    sensor_cell_map.type{idx} = 'flow';
end

% Infiltration sensors
for si = 1:sensor_config.n_infil
    idx = idx + 1;
    sensor_cell_map.row(idx) = sensor_config.infil_sensors(si, 1);
    sensor_cell_map.col(idx) = sensor_config.infil_sensors(si, 2);
    sensor_cell_map.state_idx(idx) = 6;
    sensor_cell_map.type{idx} = 'infil';
end
end

function var_y = store_sensor_variances(var_y, K_yy_local, row, col, k, sensor_config, sensor_cell_map)
% Store local K_yy variances into the global var_y array
% NOTE: Must return var_y because MATLAB passes arrays by value.

n_total = sensor_config.n_total;
local_idx = 0;

for gi = 1:n_total
    if sensor_cell_map.row(gi) == row && sensor_cell_map.col(gi) == col
        local_idx = local_idx + 1;
        if local_idx <= size(K_yy_local, 1)
            var_y(gi, k) = max(0, K_yy_local(local_idx, local_idx));
        end
    end
end
end

function y_mean = build_measurement_trajectory(h_traj, Q_l_traj, Q_r_traj, ...
    Q_u_traj, Q_d_traj, F_traj, sensor_config, sensor_cell_map, nt)
% Build the mean measurement trajectory y = C * x (no noise)

n_total = sensor_config.n_total;
y_mean = zeros(n_total, nt);

state_trajs = {h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj};

for si = 1:n_total
    r = sensor_cell_map.row(si);
    c = sensor_cell_map.col(si);
    si_idx = sensor_cell_map.state_idx(si);
    traj = state_trajs{si_idx};
    y_mean(si, :) = squeeze(traj(r, c, :))';
end
end

function R_val = get_rainfall_from_schedule(schedule, t)
% Get rainfall rate at time t from a schedule [t_start, t_end, rate]
R_val = 0;
for ri = 1:size(schedule, 1)
    if t >= schedule(ri, 1) && t < schedule(ri, 2)
        R_val = schedule(ri, 3);
        return;
    end
end
end

function merged = merge_struct_defaults(base_struct, override_struct)
% Merge override fields into base struct (defaults)
merged = base_struct;
fields = fieldnames(override_struct);
for ii = 1:numel(fields)
    merged.(fields{ii}) = override_struct.(fields{ii});
end
end

function A = build_dae_jacobian_6x6_analytical_local(y0, dt, ...
    n_cell, Ks_cell, psi_cell, eta_cell, ...
    dx, dy, h_0, h_smooth, h_flow_threshold, v_limit, F_min, has_infiltration)
% ANALYTICAL 6x6 DAE Jacobian —  approach
%
% Enhanced with regularization matching build_dae_jacobian_6x6.m:
%   - Mollified depth: h_reg = sqrt(h^2 + h_smooth^2)
%   - Flow factor smoothstep at h_flow_threshold
%   - Velocity cap detection (v_limit)
%
% ~50x faster than numerical perturbation.
%
% State order: [h, Q_l, Q_r, Q_u, Q_d, F]

cell_area = dx * dy;

h = max(y0(1), 0);
Q_l = y0(2);
Q_r = y0(3);
Q_u = y0(4);
Q_d = y0(5);
F = max(y0(6), F_min);

A = zeros(6, 6);

%% Infiltration derivatives (exact analytical)
if has_infiltration
    [~, df_dh, df_dF] = compute_infiltration_with_derivatives_local(h, F, Ks_cell, eta_cell, psi_cell, F_min);
else
    df_dh = 0;
    df_dF = 0;
end

%% Manning derivatives with mollified regularization
% Matches cell_fluxes physics in build_dae_jacobian_6x6.m

% Mollified effective depth
h_reg = sqrt(h^2 + h_smooth^2);
h_eff = max(h_reg - h_0, h_smooth);

% dh_eff/dh via chain rule through mollifier
if h_reg - h_0 >= h_smooth
    dh_eff_dh = h / h_reg;  % d(sqrt(h^2 + h_smooth^2))/dh
else
    dh_eff_dh = 0;  % clamped at h_smooth floor
end

% Flow factor: C^1 smoothstep transition at h_flow_threshold
if h > 0 && h_flow_threshold > 0
    t_flow = min(1.0, h / h_flow_threshold);
    flow_factor = t_flow^2 * (3 - 2*t_flow);
    if h < h_flow_threshold
        dflow_dh = 6 * t_flow * (1 - t_flow) / h_flow_threshold;
    else
        dflow_dh = 0;
    end
else
    flow_factor = 0;
    dflow_dh = 0;
end

% Velocity cap detection from Q values
is_capped = false;
if flow_factor > 1e-10 && h_eff > 1e-10
    v_actual = ((Q_l + Q_r)/dy + (Q_u + Q_d)/dx) / (flow_factor * h_eff);
    is_capped = (v_actual >= v_limit * 0.99);
end

% Manning power: 5/3 for power-law, 1 when velocity-limited
if is_capped
    manning_power = 1.0;
else
    manning_power = 5/3;
end

% dQ/dh = Q * (p * dh_eff_dh / h_eff + dflow_dh / flow_factor)
if flow_factor > 1e-10
    common_factor = manning_power * dh_eff_dh / h_eff + dflow_dh / flow_factor;
else
    % Near-zero flow_factor: Q ~ 0, retain h_eff derivative only
    common_factor = manning_power * dh_eff_dh / h_eff;
end

dQ_l_dh = Q_l * common_factor;
dQ_r_dh = Q_r * common_factor;
dQ_u_dh = Q_u * common_factor;
dQ_d_dh = Q_d * common_factor;

dQ_out_dh = dQ_l_dh + dQ_r_dh + dQ_u_dh + dQ_d_dh;

%% Row 1: h equation (differential - BDF)
df_h_dh = -df_dh - dQ_out_dh / cell_area;
A(1,1) = 1 - dt * df_h_dh;
A(1,2) = dt / cell_area;
A(1,3) = dt / cell_area;
A(1,4) = dt / cell_area;
A(1,5) = dt / cell_area;
A(1,6) = dt * df_dF;

%% Rows 2-5: Q algebraic constraints (NO I term)
A(2,1) = -dQ_l_dh;  A(2,2) = 1;
A(3,1) = -dQ_r_dh;  A(3,3) = 1;
A(4,1) = -dQ_u_dh;  A(4,4) = 1;
A(5,1) = -dQ_d_dh;  A(5,5) = 1;

%% Row 6: F equation (differential - BDF)
A(6,1) = -dt * df_dh;
A(6,6) = 1 - dt * df_dF;

end

function [f_infil, df_dh, df_dF] = compute_infiltration_with_derivatives_local(h, F, Ks, eta, psi_f, F_min)
% Exact derivatives for the Green-Ampt infiltration model
h_safe = max(h, 0);
F_safe = max(F, F_min);
eta_h = eta * (h_safe + psi_f);
f_capacity = Ks * (1 + eta_h / F_safe);

h_avail = 0.005;  % 5 mm transition
x_norm = min(1.0, h_safe / h_avail);
water_factor = x_norm * x_norm * (3 - 2 * x_norm);

% d(water_factor)/dh
if h_safe <= 0
    dw_dh = 0;
elseif h_safe < h_avail
    dx_dh = 1 / h_avail;
    dw_dx = 6 * x_norm * (1 - x_norm);
    dw_dh = dw_dx * dx_dh;
else
    dw_dh = 0;
end

% d(f_capacity)/dh
if h_safe <= 0
    dfc_dh = 0;
else
    dfc_dh = Ks * eta / F_safe;
end

% d(f_capacity)/dF
if F > F_min
    dfc_dF = Ks * (-eta_h / (F_safe^2));
else
    dfc_dF = 0;
end

df_dh = dfc_dh * water_factor + f_capacity * dw_dh;
df_dF = dfc_dF * water_factor;
f_infil = f_capacity * water_factor;
end
