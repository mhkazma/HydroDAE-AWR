function pse_results = run_pse_global_sparse(sim_params, sim_results, t_out, x_out, pse_config, sensor_config)
% RUN_PSE_GLOBAL_SPARSE  Global sparse PSE with measurement matrix
%
% MATHEMATICAL FRAMEWORK:
% =======================
% Unlike the cell-by-cell approach (run_pse_with_measurements), this function
% builds the FULL global sparse Jacobian A for the entire domain and solves:
%
%   K_xx = (A^T W A)^{-1}     (GLS/BLUE estimator)
%
% where A is (n_total + n_sensors) x n_total sparse, with:
%   - n_total = n_active_cells * 6 states per cell
%   - Augmented: A^s = [E; C], W = blkdiag(K_bb^{-1}, K_vv^{-1})
%
% KEY DIFFERENCE vs cell-by-cell:
%   The global Jacobian E captures CROSS-CELL COUPLING through Manning
%   routing. When cell i has a sensor, that measurement information
%   propagates to NEIGHBORING cells through the physics (shared water
%   surface slope). This gives variance reduction at UNMEASURED cells.
%
% GLOBAL JACOBIAN STRUCTURE:
%   E is block-sparse with 6x6 blocks:
%   - Diagonal blocks: dR_i/dx_i  (self-coupling, ANALYTICAL)
%   - Off-diagonal blocks: dR_i/dx_j (neighbor coupling through Manning)
%     Only entries: dQ_{dir}/dh_neighbor (rows 2-5, col 1 of neighbor block)
%
% SPARSITY:
%   For a 2D grid with 4-connectivity, each cell couples to at most 4 neighbors.
%   Total nonzeros: ~4000 * (36 + 4) = ~160,000 for a 4000-cell grid.
%   System size: 24000 x 24000 with fill ratio ~0.03%
%
% DIAGONAL EXTRACTION:
%   We don't compute the full K_xx = inv(AtWA). Instead:
%   1. Sparse Cholesky factorization of AtWA
%   2. Batch forward/backward substitution to extract diag(inv(AtWA))
%   3. For K_yy: solve AtWA \ C' (only n_sensors columns)
%
% INPUTS:
%   sim_params    - Simulation parameters
%   sim_results   - DAE simulation results (with .R rainfall)
%   t_out         - Time vector [nt x 1]
%   x_out         - State trajectory [nt x n_states_total]
%   pse_config    - PSE configuration (CVs, confidence level, etc.)
%   sensor_config - Sensor configuration from configure_sensors()
%
% OUTPUTS:
%   pse_results   - Structure matching run_pse_with_measurements output
%
% REFERENCE:
%   Theorem 1, Corollary 1
%
% Author: M. H. Kazma
% Date: February 2026

fprintf('========================================\n');
fprintf('GLOBAL SPARSE PSE\n');
fprintf('========================================\n');
fprintf('  Building FULL domain sparse Jacobian\n');
fprintf('  K_xx = (A^T W A)^{-1}  (GLS/BLUE)\n');
fprintf('  Cross-cell coupling through Manning routing\n\n');

%% ==================================================================
%  1. INITIALIZATION
%  ==================================================================

if nargin < 5 || isempty(pse_config), pse_config = struct(); end
if nargin < 6 || isempty(sensor_config)
    error('sensor_config required. Use configure_sensors() to create it.');
end

% Default PSE config
if ~isfield(pse_config, 'rainfall_cv'), pse_config.rainfall_cv = 0.20; end
if ~isfield(pse_config, 'rainfall_sigma_abs'), pse_config.rainfall_sigma_abs = 0; end
if ~isfield(pse_config, 'manning_n_cv'), pse_config.manning_n_cv = 0.20; end
if ~isfield(pse_config, 'Ks_cv'), pse_config.Ks_cv = 0.20; end
if ~isfield(pse_config, 'psi_cv'), pse_config.psi_cv = 0.20; end
if ~isfield(pse_config, 'eta_cv'), pse_config.eta_cv = 0.20; end
debug_mode = isfield(pse_config, 'debug') && pse_config.debug;

% Parameter distribution configuration (match run_pse_dae / run_pse_with_measurements)
param_dist = struct();
param_dist.rainfall = 'lognormal';
param_dist.manning_n = 'lognormal';
param_dist.Ks = 'lognormal';
param_dist.psi = 'lognormal';
param_dist.eta = 'truncated_gaussian';
if isfield(pse_config, 'param_distributions') && isstruct(pse_config.param_distributions)
    param_dist = merge_struct_defaults_local(param_dist, pse_config.param_distributions);
end

param_bounds = struct();
param_bounds.rainfall = [0, inf];
param_bounds.manning_n = [0, inf];
param_bounds.Ks = [0, inf];
param_bounds.psi = [0, inf];
param_bounds.eta = [0, 1];
if isfield(pse_config, 'param_bounds') && isstruct(pse_config.param_bounds)
    param_bounds = merge_struct_defaults_local(param_bounds, pse_config.param_bounds);
end

%% Temporal measurement schedule
% Controls WHEN sensors provide data. At measurement times, the augmented
% system A^s = [E; C] is used (over-determined → variance reduction).
% At non-measurement times, only physics A^s = E propagates uncertainty.
if isfield(pse_config, 'measurement_schedule') && isstruct(pse_config.measurement_schedule)
    meas_schedule = pse_config.measurement_schedule;
    use_temporal_schedule = true;
    fprintf('  Temporal scheduling: %s\n', meas_schedule.description);
else
    meas_schedule = [];
    use_temporal_schedule = false;
end

% Validate routing/formulation
if ~isfield(sim_params, 'routing_method') || ~strcmp(sim_params.routing_method, 'manning_ca')
    error('run_pse_global_sparse supports only manning_ca routing.');
end
if ~isfield(sim_params, 'n_states_per_cell') || sim_params.n_states_per_cell < 6
    error('run_pse_global_sparse expects 6-state formulation.');
end

%% Grid info
ny = sim_params.grid_size(1);
nx = sim_params.grid_size(2);
dx = sim_params.dx;
dy = sim_params.dy;
Z = sim_params.Z;

use_compressed = isfield(sim_params, 'use_compressed_indexing') && sim_params.use_compressed_indexing;
if use_compressed
    compressed_to_grid = sim_params.compressed_to_grid;
    grid_to_compressed = sim_params.grid_to_compressed; %#ok<NASGU>
    n_active = sim_params.valid_cells;
    active_mask = false(ny, nx);
    active_mask(sub2ind([ny, nx], compressed_to_grid(:,1), compressed_to_grid(:,2))) = true;
else
    n_active = nx * ny;
    active_mask = true(ny, nx);
    grid_to_compressed = []; %#ok<NASGU>
    compressed_to_grid = [];
end

n_states_per_cell = 6;
n_total = n_active * n_states_per_cell;
fprintf('  Domain: %d active cells, %d total states\n', n_active, n_total);

%% Physical parameters
if isfield(sim_params, 'n_manning')
    n_manning = sim_params.n_manning;
    if isscalar(n_manning), n_manning = n_manning * ones(ny, nx); end
else
    error('Manning n not found');
end

has_infiltration = isfield(sim_params, 'Ks');
if has_infiltration
    Ks = sim_params.Ks; eta = sim_params.eta;
    psi_f = sim_params.psi_f; F_min = sim_params.F_min;
else
    Ks = 0; eta = 0; psi_f = 0; F_min = 0;
end

outlet_type = 1; outlet_indices = []; slope_outlet = 0.01; outlet_direction = 'right';
if isfield(sim_params, 'outlet_type'), outlet_type = sim_params.outlet_type; end
if isfield(sim_params, 'outlet_indices'), outlet_indices = sim_params.outlet_indices; end
if isfield(sim_params, 'slope_outlet'), slope_outlet = sim_params.slope_outlet; end
if isfield(sim_params, 'outlet_direction'), outlet_direction = sim_params.outlet_direction; end

if isfield(sim_params, 'h_smooth'), h_smooth = sim_params.h_smooth;
else
    if isfield(sim_params, 'slope_max') && sim_params.slope_max > 0.2
        h_smooth = 2e-3; else, h_smooth = 1e-4; end
end
if isfield(sim_params, 'h_0'), h_0 = sim_params.h_0; else, h_0 = 0; end
if isfield(sim_params, 'h_flow_threshold'), h_flow_threshold = sim_params.h_flow_threshold; else, h_flow_threshold = 1e-3; end
if isfield(sim_params, 'v_limit'), v_limit = sim_params.v_limit; else, v_limit = 8.0; end

%% Time discretization
nt = length(t_out);

% Optional: resample trajectory to uniform time step to reduce K_xx jumps
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

%% Pre-compute variance parameters (avoid repeated computation in inner loop)

% Pre-compute sigma_ln_Ks ONCE
if has_infiltration
    if isfield(pse_config, 'sigma_ln_Ks')
        sig_lnKs = pse_config.sigma_ln_Ks;
    else
        sig_lnKs = sqrt(log(1 + pse_config.Ks_cv^2));
    end
else
    sig_lnKs = 0;
end

% PSE computation stride: compute every N timesteps
if isfield(pse_config, 'stride') && pse_config.stride > 1
    pse_stride = pse_config.stride;
else
    pse_stride = 1;
end

% Use analytical Jacobian for diagonal blocks unless user forces numerical
use_analytical_jacobian = true;
if isfield(pse_config, 'use_numerical_jacobian') && pse_config.use_numerical_jacobian
    use_analytical_jacobian = false;
end

%% ==================================================================
%  2. RAINFALL EXTRACTION (with rain_matrix support)
%  ==================================================================

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
            if use_compressed
                row_c = compressed_to_grid(cc, 1);
                col_c = compressed_to_grid(cc, 2);
            else
                row_c = ceil(cc / nx);
                col_c = mod(cc - 1, nx) + 1;
            end
            R_rain_matrix_grid(row_c, col_c, kk) = rm(cc, t_idx);
        end
    end
    R_time = squeeze(mean(mean(R_rain_matrix_grid, 1), 2))';
    fprintf('  PSE: rain_matrix mode (%d x %d x %d)\n', ny, nx, nt);
elseif isstruct(R_in)
    if isfield(R_in, 'schedule')
        schedule = R_in.schedule;
        R_time = zeros(1, nt);
        for kk = 1:nt
            R_time(kk) = get_rain_schedule_local(schedule, t_out(kk));
        end
        if isfield(R_in, 'pattern') && ~isempty(R_in.pattern)
            [R_pattern_scalar, R_pattern_grid] = build_rain_pattern_local(R_in.pattern, ...
                use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
        end
    elseif isfield(R_in, 'base') && isfield(R_in, 'duration')
        base = R_in.base;
        [R_pattern_scalar, R_pattern_grid] = build_rain_pattern_local(base, ...
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
    [R_pattern_scalar, R_pattern_grid] = build_rain_pattern_local(R_in, ...
        use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
    R_time = ones(1, nt);
end

if isempty(R_time)
    R_time = ones(1, nt);
end

%% ==================================================================
%  3. EXTRACT DETERMINISTIC TRAJECTORIES
%  ==================================================================

fprintf('  Extracting trajectories...\n');
[h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj] = ...
    extract_trajectories_global(x_out, nt, ny, nx, n_states_per_cell, ...
    use_compressed, n_active, compressed_to_grid);

%% ==================================================================
%  4. BUILD CELL INDEXING AND NEIGHBOR MAP
%  ==================================================================

cell_rows = zeros(n_active, 1);
cell_cols = zeros(n_active, 1);
cell_index_map = zeros(ny, nx);
if use_compressed
    cell_rows = compressed_to_grid(:, 1);
    cell_cols = compressed_to_grid(:, 2);
    for a = 1:n_active
        cell_index_map(cell_rows(a), cell_cols(a)) = a;
    end
else
    a = 0;
    for i = 1:ny
        for j = 1:nx
            a = a + 1;
            cell_rows(a) = i;
            cell_cols(a) = j;
            cell_index_map(i, j) = a;
        end
    end
end

% Neighbor map: for each cell, indices of [left, right, up, down] neighbors
% 0 means no neighbor (boundary or inactive)
neighbor_ci = zeros(n_active, 4);
for a = 1:n_active
    r = cell_rows(a);
    c = cell_cols(a);
    if c > 1  && active_mask(r, c-1), neighbor_ci(a, 1) = cell_index_map(r, c-1); end
    if c < nx && active_mask(r, c+1), neighbor_ci(a, 2) = cell_index_map(r, c+1); end
    if r > 1  && active_mask(r-1, c), neighbor_ci(a, 3) = cell_index_map(r-1, c); end
    if r < ny && active_mask(r+1, c), neighbor_ci(a, 4) = cell_index_map(r+1, c); end
end

fprintf('  Neighbor map built: avg %.1f neighbors/cell\n', ...
    mean(sum(neighbor_ci > 0, 2)));

%% ==================================================================
%  5. BUILD GLOBAL MEASUREMENT MATRIX C (time-independent, sparse)
%  ==================================================================

n_sensors = sensor_config.n_total;
[C_global, K_vv_vec] = build_global_C_matrix_local(sensor_config, cell_index_map, ...
    n_total, n_sensors, n_states_per_cell);

K_vv_global = spdiags(K_vv_vec, 0, n_sensors, n_sensors);

fprintf('  Global C matrix: %d x %d, nnz = %d\n', size(C_global,1), size(C_global,2), nnz(C_global));

%% ==================================================================
%  6. PREALLOCATE OUTPUTS
%  ==================================================================

var_h = zeros(ny, nx, nt);
var_F = zeros(ny, nx, nt);
var_Ql = zeros(ny, nx, nt);
var_Qr = zeros(ny, nx, nt);
var_Qu = zeros(ny, nx, nt);
var_Qd = zeros(ny, nx, nt);
var_y = zeros(n_sensors, nt);

%% ==================================================================
%  7. MAIN TIME LOOP (stride-based)
%  ==================================================================

if use_analytical_jacobian
    fprintf('  Computing K_xx with ANALYTICAL diagonal Jacobian + numerical cross-cell...\n');
else
    fprintf('  Computing K_xx with NUMERICAL Jacobian (full perturbation)...\n');
end
fprintf('  Sensors: %d depth, %d flow, %d infil (total %d)\n', ...
    sensor_config.n_depth, sensor_config.n_flow, sensor_config.n_infil, n_sensors);
if use_temporal_schedule
    fprintf('  Temporal schedule: %d/%d timesteps with measurements (%.1f%%)\n', ...
        meas_schedule.n_measurements, meas_schedule.n_total, 100*meas_schedule.fraction);
else
    fprintf('  Temporal schedule: continuous (measurements at every timestep)\n');
end
fprintf('  Stride=%d, %d active cells, %d timesteps\n', pse_stride, n_active, nt);

tic;
wall_start = tic;

% Batch size for diagonal extraction
batch_size = min(500, n_total);

% Determine which timesteps to compute PSE for
pse_timesteps = 2:pse_stride:nt;
if pse_timesteps(end) ~= nt
    pse_timesteps = [pse_timesteps, nt];
end
n_pse_steps = length(pse_timesteps);
pse_step_count = 0;

% Variance reduction tracking (for debug diagnostics)
variance_reduction_sum = 0;
variance_reduction_count = 0;

for ki = 1:n_pse_steps
    k = pse_timesteps(ki);

    % Previous PSE timestep (for variance retrieval)
    if ki > 1
        k_prev_pse = pse_timesteps(ki-1);
    else
        k_prev_pse = 1;
    end

    % For stride > 1, accumulate dt across skipped steps
    if pse_stride > 1 && ki > 1
        dt = sum(dt_vec(k_prev_pse:k-1));
    else
        dt = dt_vec(k-1);
    end

    % ETA progress reporting
    pse_step_count = pse_step_count + 1;
    if mod(pse_step_count, 10) == 0 || pse_step_count == 1 || k == nt
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

    %% 7a. Build global sparse Jacobian A (n_total x n_total)
    A_global = build_global_jacobian_at_timestep_local(k, dt, ...
        h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj, ...
        Z, n_manning, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, ...
        active_mask, has_infiltration, Ks, eta, psi_f, F_min, ...
        outlet_type, outlet_indices, outlet_direction, slope_outlet, ...
        n_active, n_total, n_states_per_cell, neighbor_ci, ...
        cell_rows, cell_cols, ny, nx, use_analytical_jacobian);

    %% 7b. Build global block-diagonal K_bb (n_total x n_total sparse)
    K_bb_global = build_global_kbb_at_timestep_local(k, dt, k_prev_pse, ...
        h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, ...
        var_h, var_F, ...
        R_time, R_pattern_scalar, R_pattern_grid, ...
        use_rain_matrix_pse, R_rain_matrix_grid, ...
        n_manning, has_infiltration, Ks, eta, psi_f, F_min, sig_lnKs, ...
        n_active, n_total, n_states_per_cell, ...
        cell_rows, cell_cols, ny, nx, pse_config);

    %% 7c. Augment system: A^s = [A; C], K_bb^s = blkdiag(K_bb, K_vv)
    % Temporal scheduling: only include measurements at scheduled times
    if use_temporal_schedule
        use_meas_this_step = meas_schedule.is_measurement(k) && n_sensors > 0;
    else
        use_meas_this_step = n_sensors > 0;
    end

    if use_meas_this_step
        A_aug = [A_global; C_global];
    else
        A_aug = A_global;
    end

    %% 7d. Build W = inv(K_bb_aug) — block-diagonal, sparse
    if use_meas_this_step
        W_global = build_global_W_inverse_local(K_bb_global, K_vv_vec, ...
            n_active, n_sensors, n_total, n_states_per_cell);
    else
        % No measurements: W = inv(K_bb) only (square system)
        W_global = build_global_W_inverse_local(K_bb_global, [], ...
            n_active, 0, n_total, n_states_per_cell);
    end

    %% 7e. Fisher information: AtWA = A^T W A
    WA = W_global * A_aug;
    AtWA = A_aug' * WA;
    AtWA = 0.5 * (AtWA + AtWA');  % Enforce symmetry

    %% 7f. Extract diagonal of K_xx = inv(AtWA)
    diag_Kxx = extract_diagonal_of_sparse_inverse_local(AtWA, batch_size);

    %% 7g. Compute K_yy for sensors (only at measurement times)
    if use_meas_this_step
        try
            X_sensors = AtWA \ C_global';
            K_yy_full = C_global * X_sensors + K_vv_global;
            var_y(:, k) = max(0, full(diag(K_yy_full)));
        catch
            lambda = 1e-10 * max(abs(diag(AtWA)));
            AtWA_reg = AtWA + lambda * speye(n_total);
            X_sensors = AtWA_reg \ C_global';
            K_yy_full = C_global * X_sensors + K_vv_global;
            var_y(:, k) = max(0, full(diag(K_yy_full)));
        end
    end

    %% 7h. Store variances per cell
    for a = 1:n_active
        r = cell_rows(a);
        c = cell_cols(a);
        base = (a - 1) * n_states_per_cell;
        var_h(r, c, k)  = max(0, diag_Kxx(base + 1));
        var_Ql(r, c, k) = max(0, diag_Kxx(base + 2));
        var_Qr(r, c, k) = max(0, diag_Kxx(base + 3));
        var_Qu(r, c, k) = max(0, diag_Kxx(base + 4));
        var_Qd(r, c, k) = max(0, diag_Kxx(base + 5));
        var_F(r, c, k)  = max(0, diag_Kxx(base + 6));
    end

    %% 7i. Variance reduction tracking (debug mode, measurement timesteps)
    if debug_mode && use_meas_this_step
        % Compute prior diag (without measurements) for comparison
        W_prior = build_global_W_inverse_local(K_bb_global, [], ...
            n_active, 0, n_total, n_states_per_cell);
        WA_prior = W_prior * A_global;
        AtWA_prior = A_global' * WA_prior;
        AtWA_prior = 0.5 * (AtWA_prior + AtWA_prior');
        diag_Kxx_prior = extract_diagonal_of_sparse_inverse_local(AtWA_prior, batch_size);

        % Compare posterior/prior at sensor locations (h state only)
        for si = 1:n_sensors
            s_row = sensor_config.rows(si);
            s_col = sensor_config.cols(si);
            s_ci = cell_index_map(s_row, s_col);
            if s_ci > 0
                base_s = (s_ci - 1) * n_states_per_cell;
                var_post = max(0, diag_Kxx(base_s + 1));       % h posterior
                var_prior = max(0, diag_Kxx_prior(base_s + 1)); % h prior
                if var_prior > 1e-20
                    reduction = var_post / var_prior;
                    variance_reduction_sum = variance_reduction_sum + reduction;
                    variance_reduction_count = variance_reduction_count + 1;
                end
            end
        end
    end
end

%% ==================================================================
%  8. STRIDE INTERPOLATION
%  ==================================================================

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
                    var_h(:,:,kk)  = var_h(:,:,k_prev);
                    var_F(:,:,kk)  = var_F(:,:,k_prev);
                    var_Ql(:,:,kk) = var_Ql(:,:,k_prev);
                    var_Qr(:,:,kk) = var_Qr(:,:,k_prev);
                    var_Qu(:,:,kk) = var_Qu(:,:,k_prev);
                    var_Qd(:,:,kk) = var_Qd(:,:,k_prev);
                    if n_sensors > 0
                        var_y(:,kk) = var_y(:,k_prev);
                    end
                end
            end
        end
    elseif strcmp(stride_interp, 'pchip')
        % PCHIP interpolation per-cell (shape-preserving, C^1)
        t_pse = t_out(pse_timesteps);
        for ci = 1:n_active
            r = cell_rows(ci); c = cell_cols(ci);
            vals_h  = squeeze(var_h(r,c,pse_timesteps));
            vals_F  = squeeze(var_F(r,c,pse_timesteps));
            vals_Ql = squeeze(var_Ql(r,c,pse_timesteps));
            vals_Qr = squeeze(var_Qr(r,c,pse_timesteps));
            vals_Qu = squeeze(var_Qu(r,c,pse_timesteps));
            vals_Qd = squeeze(var_Qd(r,c,pse_timesteps));
            for ki = 2:n_pse_steps
                k_curr = pse_timesteps(ki);
                k_prev = pse_timesteps(ki-1);
                if k_curr - k_prev > 1
                    kk_range = (k_prev+1):(k_curr-1);
                    t_interp = t_out(kk_range);
                    var_h(r,c,kk_range)  = max(0, interp1(t_pse, vals_h, t_interp, 'pchip'));
                    var_F(r,c,kk_range)  = max(0, interp1(t_pse, vals_F, t_interp, 'pchip'));
                    var_Ql(r,c,kk_range) = max(0, interp1(t_pse, vals_Ql, t_interp, 'pchip'));
                    var_Qr(r,c,kk_range) = max(0, interp1(t_pse, vals_Qr, t_interp, 'pchip'));
                    var_Qu(r,c,kk_range) = max(0, interp1(t_pse, vals_Qu, t_interp, 'pchip'));
                    var_Qd(r,c,kk_range) = max(0, interp1(t_pse, vals_Qd, t_interp, 'pchip'));
                end
            end
        end
        % Interpolate measurement variance separately
        if n_sensors > 0
            for si = 1:n_sensors
                vals_y = var_y(si, pse_timesteps);
                for ki = 2:n_pse_steps
                    k_curr = pse_timesteps(ki);
                    k_prev = pse_timesteps(ki-1);
                    if k_curr - k_prev > 1
                        kk_range = (k_prev+1):(k_curr-1);
                        t_interp = t_out(kk_range);
                        var_y(si, kk_range) = max(0, interp1(t_pse, vals_y, t_interp, 'pchip'));
                    end
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
                    var_h(:,:,kk)  = (1-alpha) * var_h(:,:,k_prev)  + alpha * var_h(:,:,k_curr);
                    var_F(:,:,kk)  = (1-alpha) * var_F(:,:,k_prev)  + alpha * var_F(:,:,k_curr);
                    var_Ql(:,:,kk) = (1-alpha) * var_Ql(:,:,k_prev) + alpha * var_Ql(:,:,k_curr);
                    var_Qr(:,:,kk) = (1-alpha) * var_Qr(:,:,k_prev) + alpha * var_Qr(:,:,k_curr);
                    var_Qu(:,:,kk) = (1-alpha) * var_Qu(:,:,k_prev) + alpha * var_Qu(:,:,k_curr);
                    var_Qd(:,:,kk) = (1-alpha) * var_Qd(:,:,k_prev) + alpha * var_Qd(:,:,k_curr);
                    if n_sensors > 0
                        var_y(:,kk) = (1-alpha) * var_y(:,k_prev) + alpha * var_y(:,k_curr);
                    end
                end
            end
        end
    end
end

pse_wall_time = toc(wall_start);
computation_time = toc;

%% ==================================================================
%  9. SUMMARY AND DIAGNOSTICS
%  ==================================================================

fprintf('  Global sparse PSE complete in %.1f seconds (%.1f s/step)\n', ...
    pse_wall_time, pse_wall_time / max(n_pse_steps, 1));
if use_analytical_jacobian
    fprintf('  - Method: analytical diagonal + numerical cross-cell, stride=%d\n', pse_stride);
else
    fprintf('  - Method: numerical Jacobian (full), stride=%d\n', pse_stride);
end
fprintf('  - h, F: Temporal propagation + parameter injection\n');
fprintf('  - Q: Parameter injection only (algebraic, no temporal memory)\n');

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

%% ==================================================================
%  10. PACKAGE RESULTS
%  ==================================================================

h_std = sqrt(var_h);
F_std = sqrt(var_F);
Q_l_std = sqrt(var_Ql);
Q_r_std = sqrt(var_Qr);
Q_u_std = sqrt(var_Qu);
Q_d_std = sqrt(var_Qd);

pse_results = struct();

% Standard deviations
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

% Full variance arrays
pse_results.var_h = var_h;
pse_results.var_F = var_F;
pse_results.var_Ql = var_Ql;
pse_results.var_Qr = var_Qr;
pse_results.var_Qu = var_Qu;
pse_results.var_Qd = var_Qd;

% Measurement output uncertainty
pse_results.var_y = var_y;
pse_results.y_std = sqrt(max(0, var_y));
pse_results.sensor_config = sensor_config;

%% Confidence intervals
dist_type = 'gaussian';
conf_level = 0.95;
if isfield(pse_config, 'distribution_type'), dist_type = pse_config.distribution_type; end
if isfield(pse_config, 'confidence_level'), conf_level = pse_config.confidence_level; end

fprintf('  Computing %s confidence intervals (%.0f%%)...\n', dist_type, conf_level*100);

% Optional CI parameters (match run_pse_dae / run_pse_with_measurements)
ci_opts = {};
if isfield(pse_config, 'min_depth_for_gaussian')
    ci_opts = [ci_opts, {'MinDepth', pse_config.min_depth_for_gaussian}];
end
if isfield(pse_config, 'truncation_lower')
    ci_opts = [ci_opts, {'TruncLower', pse_config.truncation_lower}];
end

[pse_results.h_lower, pse_results.h_upper] = compute_confidence_intervals(...
    h_traj, h_std, dist_type, conf_level, ci_opts{:});
pse_results.h_lower = max(pse_results.h_lower, 0);  % h >= 0 (physical)

[pse_results.F_lower, pse_results.F_upper] = compute_confidence_intervals(...
    F_traj, F_std, dist_type, conf_level, ci_opts{:});
pse_results.F_lower = max(pse_results.F_lower, 0);  % F >= 0 (physical)

% Q's use Gaussian (can be negative)
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

% Measurement output CIs
if n_sensors > 0
    y_mean = build_measurement_trajectory_global_local(h_traj, Q_l_traj, Q_r_traj, ...
        Q_u_traj, Q_d_traj, F_traj, sensor_config, cell_index_map, nt);
    pse_results.y_mean = y_mean;
    if strcmp(dist_type, 'gaussian')
        z_alpha = norminv(0.5 + conf_level/2);
    else
        z_alpha = 1.96;
    end
    pse_results.y_lower = y_mean - z_alpha * pse_results.y_std;
    pse_results.y_upper = y_mean + z_alpha * pse_results.y_std;
end

%% Metadata (match run_pse_with_measurements output)
pse_results.distribution_type = dist_type;
pse_results.confidence_level = conf_level;
pse_results.param_distributions = param_dist;
pse_results.param_bounds = param_bounds;
pse_results.t = t_out;
pse_results.computation_time = computation_time;
pse_results.method = sprintf('PSE GLOBAL SPARSE (cross-cell coupling) - %s', dist_type);
pse_results.n_cells = n_active;
pse_results.n_timesteps = nt;
pse_results.grid_size = [ny, nx];
pse_results.n_sensors = n_sensors;
pse_results.n_total_states = n_total;

%% Variance reduction analysis (matches run_pse_with_measurements output)
if debug_mode && variance_reduction_count > 0
    avg_reduction = variance_reduction_sum / variance_reduction_count;
    fprintf('  Average variance reduction at measured cells: %.1f%%\n', ...
        (1 - avg_reduction) * 100);
    pse_results.avg_variance_reduction = avg_reduction;
end

fprintf('  PSE completed in %.2f seconds\n', computation_time);
fprintf('========================================\n');

end % main function

%% ========================================================================
%  GLOBAL JACOBIAN ASSEMBLY
%  ========================================================================

function A = build_global_jacobian_at_timestep_local(k, dt, ...
    h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj, ...
    Z, n_manning, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, ...
    active_mask, has_infiltration, Ks, eta, psi_f, F_min, ...
    outlet_type, outlet_indices, outlet_direction, slope_outlet, ...
    n_active, n_total, nsp, neighbor_ci, ...
    cell_rows, cell_cols, ny, nx, use_analytical_jacobian) %#ok<INUSL>
% Build the global sparse Jacobian at timestep k
%
% The Jacobian has:
%   - Diagonal 6×6 blocks: dR_i/dx_i (self-coupling)
%     → ANALYTICAL when use_analytical_jacobian = true (~50x faster)
%   - Off-diagonal entries: dR_i/dh_j (neighbor coupling through Manning)
%     → Always NUMERICAL (WSE slope depends on neighbor h)

% Preallocate triplets
max_nnz = n_active * 36 + n_active * 4 * 4;  % diag + off-diag
I_idx = zeros(max_nnz, 1);
J_idx = zeros(max_nnz, 1);
V_val = zeros(max_nnz, 1);
ptr = 0;

for a = 1:n_active
    row = cell_rows(a);
    col = cell_cols(a);
    base = (a - 1) * nsp;

    % Operating point
    h_curr = h_traj(row, col, k);
    F_curr = F_traj(row, col, k);
    Q_l = Q_l_traj(row, col, k);
    Q_r = Q_r_traj(row, col, k);
    Q_u = Q_u_traj(row, col, k);
    Q_d = Q_d_traj(row, col, k);
    y0 = [h_curr; Q_l; Q_r; Q_u; Q_d; F_curr];

    n_cell = n_manning(row, col);

    if has_infiltration
        Ks_cell = get_param_local(Ks, row, col);
        psi_cell = get_param_local(psi_f, row, col);
        eta_cell = get_param_local(eta, row, col);
    else
        Ks_cell = 0; psi_cell = 0; eta_cell = 0;
    end

    %% Diagonal block: local 6×6 Jacobian (self-coupling)
    if use_analytical_jacobian
        % Analytical: regularized dQ/dh with mollifier, flow_factor, v_limit
        A_local = build_dae_jacobian_6x6_analytical_local(y0, dt, ...
            n_cell, Ks_cell, psi_cell, eta_cell, ...
            dx, dy, h_0, h_smooth, h_flow_threshold, v_limit, F_min, has_infiltration);
    else
        % Numerical: central difference perturbation
        A_local = build_dae_jacobian_6x6_numerical_local(y0, row, col, k, dt, ...
            n_cell, Ks_cell, psi_cell, eta_cell, ...
            h_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, ...
            active_mask, F_min, has_infiltration, ...
            outlet_type, outlet_indices, outlet_direction, slope_outlet);
    end

    % Store diagonal block entries
    for ii = 1:nsp
        for jj = 1:nsp
            if abs(A_local(ii, jj)) > 1e-30
                ptr = ptr + 1;
                I_idx(ptr) = base + ii;
                J_idx(ptr) = base + jj;
                V_val(ptr) = A_local(ii, jj);
            end
        end
    end

    %% Off-diagonal: neighbor coupling dQ_dir/dh_neighbor
    % Only Q equations (rows 2-5) depend on neighbor h
    if h_curr < 1e-10
        continue;  % No flow at dry cells — skip neighbor coupling
    end

    h_eps = max(1e-6, 1e-3 * max(h_curr, 1e-6));
    h_grid_k = h_traj(:, :, k);

    for dir = 1:4
        nci = neighbor_ci(a, dir);
        if nci == 0, continue; end

        nr = cell_rows(nci);
        nc = cell_cols(nci);
        neighbor_base = (nci - 1) * nsp;
        h_neighbor = h_traj(nr, nc, k);

        % Perturb neighbor h and compute flux change at THIS cell
        h_grid_plus = h_grid_k;
        h_grid_plus(nr, nc) = h_neighbor + h_eps;
        [Ql_p, Qr_p, Qu_p, Qd_p] = cell_fluxes_global_local(h_grid_plus, row, col, ...
            n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask);

        h_grid_minus = h_grid_k;
        h_grid_minus(nr, nc) = max(h_neighbor - h_eps, 0);
        [Ql_m, Qr_m, Qu_m, Qd_m] = cell_fluxes_global_local(h_grid_minus, row, col, ...
            n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask);

        eps_actual = h_grid_plus(nr, nc) - h_grid_minus(nr, nc);
        if eps_actual < 1e-15, continue; end

        % dQ/dh_neighbor for each Q equation row
        dQl_dhn = -(Ql_p - Ql_m) / eps_actual;
        dQr_dhn = -(Qr_p - Qr_m) / eps_actual;
        dQu_dhn = -(Qu_p - Qu_m) / eps_actual;
        dQd_dhn = -(Qd_p - Qd_m) / eps_actual;

        deriv_vals = [dQl_dhn, dQr_dhn, dQu_dhn, dQd_dhn];
        q_rows = [2, 3, 4, 5];

        for qi = 1:4
            if abs(deriv_vals(qi)) > 1e-30
                ptr = ptr + 1;
                if ptr > max_nnz
                    max_nnz = max_nnz * 2;
                    I_idx(end+1:max_nnz) = 0;
                    J_idx(end+1:max_nnz) = 0;
                    V_val(end+1:max_nnz) = 0;
                end
                I_idx(ptr) = base + q_rows(qi);
                J_idx(ptr) = neighbor_base + 1;  % h of neighbor
                V_val(ptr) = deriv_vals(qi);
            end
        end
    end
end

I_idx = I_idx(1:ptr);
J_idx = J_idx(1:ptr);
V_val = V_val(1:ptr);
A = sparse(I_idx, J_idx, V_val, n_total, n_total);
end

%% ========================================================================
%  ANALYTICAL 6×6 DAE JACOBIAN (diagonal blocks)
%  ========================================================================

function A = build_dae_jacobian_6x6_analytical_local(y0, dt, ...
    n_cell, Ks_cell, psi_cell, eta_cell, ...
    dx, dy, h_0, h_smooth, h_flow_threshold, v_limit, F_min, has_infiltration) %#ok<INUSD>
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

%% Rows 2-5: Q algebraic constraints (NO BDF I term!)
A(2,1) = -dQ_l_dh;  A(2,2) = 1;
A(3,1) = -dQ_r_dh;  A(3,3) = 1;
A(4,1) = -dQ_u_dh;  A(4,4) = 1;
A(5,1) = -dQ_d_dh;  A(5,5) = 1;

%% Row 6: F equation (differential - BDF)
A(6,1) = -dt * df_dh;
A(6,6) = 1 - dt * df_dF;

end

%% ========================================================================
%  NUMERICAL 6×6 DAE JACOBIAN (diagonal blocks — fallback)
%  ========================================================================

function A = build_dae_jacobian_6x6_numerical_local(y0, row, col, k, dt, ...
    n_cell, Ks_cell, psi_cell, eta_cell, ...
    h_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, ...
    active_mask, F_min, has_infiltration, ...
    outlet_type, outlet_indices, outlet_direction, slope_outlet) %#ok<INUSL>
% NUMERICAL 6×6 DAE Jacobian via central difference perturbation
% Used as fallback when use_analytical_jacobian = false

cell_area = dx * dy;
h = max(y0(1), 0);
F = max(y0(6), F_min);

A = zeros(6, 6);

%% Infiltration derivatives (exact)
if has_infiltration
    [~, df_dh, df_dF] = compute_infiltration_with_derivatives_local(h, F, Ks_cell, eta_cell, psi_cell, F_min);
else
    df_dh = 0;
    df_dF = 0;
end

%% Manning derivatives (numerical)
h_eps = max(1e-6, 1e-3 * max(h, 1e-6));
h_grid = h_traj(:,:,k);

h_grid(row, col) = h + h_eps;
[Qlp, Qrp, Qup, Qdp] = cell_fluxes_global_local(h_grid, row, col, ...
    n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask);

h_grid(row, col) = max(h - h_eps, 0);
[Qlm, Qrm, Qum, Qdm] = cell_fluxes_global_local(h_grid, row, col, ...
    n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask);

% Apply outlet boundary override
is_outlet = false;
if ~isempty(outlet_indices)
    for oi = 1:size(outlet_indices, 1)
        if outlet_indices(oi,1) == row && outlet_indices(oi,2) == col
            is_outlet = true;
            break;
        end
    end
end

if is_outlet && outlet_type == 1
    h_eff_out = max(h - h_0, 1e-6);
    if strcmp(outlet_direction, 'right')
        Qrp = (1/n_cell) * (h_eff_out + h_eps)^(5/3) * sqrt(slope_outlet) * dy;
        Qrm = (1/n_cell) * max(h_eff_out - h_eps, 1e-6)^(5/3) * sqrt(slope_outlet) * dy;
    end
end

dQ_l_dh = (Qlp - Qlm) / (2 * h_eps);
dQ_r_dh = (Qrp - Qrm) / (2 * h_eps);
dQ_u_dh = (Qup - Qum) / (2 * h_eps);
dQ_d_dh = (Qdp - Qdm) / (2 * h_eps);

dQ_out_dh = dQ_l_dh + dQ_r_dh + dQ_u_dh + dQ_d_dh;

%% Row 1: h equation (BDF)
df_h_dh = -df_dh - dQ_out_dh / cell_area;
A(1,1) = 1 - dt * df_h_dh;
A(1,2) = dt / cell_area;
A(1,3) = dt / cell_area;
A(1,4) = dt / cell_area;
A(1,5) = dt / cell_area;
A(1,6) = -dt * (-df_dF);

%% Rows 2-5: Q algebraic (NO I term)
A(2,1) = -dQ_l_dh;  A(2,2) = 1;
A(3,1) = -dQ_r_dh;  A(3,3) = 1;
A(4,1) = -dQ_u_dh;  A(4,4) = 1;
A(5,1) = -dQ_d_dh;  A(5,5) = 1;

%% Row 6: F equation (BDF)
A(6,1) = -dt * df_dh;
A(6,6) = 1 - dt * df_dF;

end

%% ========================================================================
%  INFILTRATION WITH EXACT DERIVATIVES
%  ========================================================================

function [f_infil, df_dh, df_dF] = compute_infiltration_with_derivatives_local(h, F, Ks, eta, psi_f, F_min)
% Exact derivatives for Green-Ampt infiltration with smoothed water_factor
h_safe = max(h, 0);
F_safe = max(F, F_min);
eta_h = eta * (h_safe + psi_f);
f_capacity = Ks * (1 + eta_h / F_safe);

% Smoothed Hermite transition at h_avail = 5 mm
h_avail = 0.005;
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

%% ========================================================================
%  GLOBAL K_bb ASSEMBLY (block-diagonal)
%  ========================================================================

function K_bb = build_global_kbb_at_timestep_local(k, dt, k_prev_pse, ...
    h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, ...
    var_h, var_F, ...
    R_time, R_pattern_scalar, R_pattern_grid, ...
    use_rain_matrix_pse, R_rain_matrix_grid, ...
    n_manning, has_infiltration, Ks, eta, psi_f, F_min, sig_lnKs, ...
    n_active, n_total, nsp, ...
    cell_rows, cell_cols, ny, nx, pse_config) %#ok<INUSL>
% Build global block-diagonal K_bb with all improvements:
%   - rain_matrix support
%   - pre-computed sig_lnKs
%   - proper k_prev_pse for stride support
%   - dry cell fast path

max_nnz = n_active * 36;
I_idx = zeros(max_nnz, 1);
J_idx = zeros(max_nnz, 1);
V_val = zeros(max_nnz, 1);
ptr = 0;

for a = 1:n_active
    row = cell_rows(a);
    col = cell_cols(a);
    base = (a - 1) * nsp;

    h_curr = h_traj(row, col, k);
    F_curr = F_traj(row, col, k);
    Q_l = Q_l_traj(row, col, k);
    Q_r = Q_r_traj(row, col, k);
    Q_u = Q_u_traj(row, col, k);
    Q_d = Q_d_traj(row, col, k);
    n_cell = n_manning(row, col);

    % Rainfall at this cell and timestep (with rain_matrix support)
    if use_rain_matrix_pse
        R_cell = R_rain_matrix_grid(row, col, k);
    elseif ~isempty(R_pattern_grid)
        R_cell = R_pattern_grid(row, col) * R_time(k);
    else
        R_cell = R_pattern_scalar * R_time(k);
    end

    % Previous-step variances (using k_prev_pse for stride support)
    var_h_prev = var_h(row, col, k_prev_pse);
    var_F_prev = var_F(row, col, k_prev_pse);

    %% Dry cell fast path
    if h_curr < 1e-8
        K_bb_local = zeros(6, 6);

        % Rainfall uncertainty even at dry cells
        var_R_dry = (dt * pse_config.rainfall_cv * R_cell)^2 + ...
                    (dt * pse_config.rainfall_sigma_abs)^2;

        % h: small from rainfall or decaying
        if R_cell > 0
            if has_infiltration
                Ks_cell = get_param_local(Ks, row, col);
                var_Ks_dry = (Ks_cell * sig_lnKs * dt)^2;
                K_bb_local(1,1) = max(var_h_prev * 0.5, var_R_dry + var_Ks_dry);
                K_bb_local(6,6) = var_F_prev + var_R_dry + var_Ks_dry;
            else
                K_bb_local(1,1) = var_R_dry;
                K_bb_local(6,6) = var_F_prev + var_R_dry;
            end
        else
            K_bb_local(1,1) = var_h_prev * 0.1;
            if has_infiltration
                Ks_cell = get_param_local(Ks, row, col);
                K_bb_local(6,6) = var_F_prev + (Ks_cell * sig_lnKs * dt)^2;
            else
                K_bb_local(6,6) = var_F_prev;
            end
        end
        % Q: zero variance for dry cells (Q ≈ 0)

        % Store sparse entries
        for ii = 1:nsp
            for jj = 1:nsp
                if abs(K_bb_local(ii, jj)) > 1e-30
                    ptr = ptr + 1;
                    I_idx(ptr) = base + ii;
                    J_idx(ptr) = base + jj;
                    V_val(ptr) = K_bb_local(ii, jj);
                end
            end
        end
        continue;
    end

    %% Wet cell path — full K_bb construction
    if has_infiltration
        Ks_cell = get_param_local(Ks, row, col);
        psi_cell = get_param_local(psi_f, row, col);
        eta_cell = get_param_local(eta, row, col);
    else
        Ks_cell = 0; psi_cell = 0; eta_cell = 0;
    end

    K_bb_local = build_kbb_local_func(h_curr, F_curr, Q_l, Q_r, Q_u, Q_d, ...
        var_h_prev, var_F_prev, dt, R_cell, n_cell, ...
        Ks_cell, psi_cell, eta_cell, F_min, has_infiltration, pse_config, sig_lnKs);

    % Store in triplet format
    for ii = 1:nsp
        for jj = 1:nsp
            if abs(K_bb_local(ii, jj)) > 1e-30
                ptr = ptr + 1;
                I_idx(ptr) = base + ii;
                J_idx(ptr) = base + jj;
                V_val(ptr) = K_bb_local(ii, jj);
            end
        end
    end
end

I_idx = I_idx(1:ptr);
J_idx = J_idx(1:ptr);
V_val = V_val(1:ptr);
K_bb = sparse(I_idx, J_idx, V_val, n_total, n_total);
end

%% ========================================================================
%  LOCAL K_bb BUILDER (per-cell, with pre-computed sig_lnKs)
%  ========================================================================

function K_bb = build_kbb_local_func(h_curr, F_curr, Q_l, Q_r, Q_u, Q_d, ...
    var_h_prev, var_F_prev, dt, R_cell, n_cell, ...
    Ks_cell, psi_cell, eta_cell, F_min, has_infiltration, pse_config, sig_lnKs)
% Build local 6×6 K_bb with pre-computed sig_lnKs (matches run_pse_with_measurements)

K_bb = zeros(6, 6);

% Temporal propagation: C * K_xx^n * C' where C = diag([1,0,0,0,0,1])
K_bb(1,1) = var_h_prev;
K_bb(6,6) = var_F_prev;

% Rainfall uncertainty (ALWAYS include absolute floor even when R=0)
var_R = (dt * pse_config.rainfall_cv * R_cell)^2 + (dt * pse_config.rainfall_sigma_abs)^2;
K_bb(1,1) = K_bb(1,1) + var_R;

% Infiltration uncertainty
if has_infiltration && F_curr > F_min
    f_infil = Ks_cell * (1 + eta_cell * (h_curr + psi_cell) / max(F_curr, F_min));

    % Use pre-computed sig_lnKs (log-transform for Ks)
    df_dlnKs = f_infil;
    var_f_from_Ks = df_dlnKs^2 * sig_lnKs^2;

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
K_bb(1,6) = -var_f_infil;
K_bb(6,1) = -var_f_infil;

% Manning n uncertainty → creates Q-Q correlations
cv_n = pse_config.manning_n_cv;
var_n = (cv_n * n_cell)^2;
dQl_dn = -Q_l / max(abs(n_cell), 1e-10);
dQr_dn = -Q_r / max(abs(n_cell), 1e-10);
dQu_dn = -Q_u / max(abs(n_cell), 1e-10);
dQd_dn = -Q_d / max(abs(n_cell), 1e-10);

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

%% ========================================================================
%  W = inv(K_bb_aug) — BLOCK-DIAGONAL INVERSE
%  ========================================================================

function W = build_global_W_inverse_local(K_bb_global, K_vv_vec, ...
    n_active, n_sensors, n_total, nsp)
% Build W = inv(blkdiag(K_bb, K_vv)) as a sparse matrix

max_nnz = n_active * 36 + n_sensors;
I_idx = zeros(max_nnz, 1);
J_idx = zeros(max_nnz, 1);
V_val = zeros(max_nnz, 1);
ptr = 0;

for a = 1:n_active
    base = (a - 1) * nsp;
    idx = base + (1:nsp);

    K_block = full(K_bb_global(idx, idx));

    % Regularize if needed
    rc = rcond(K_block);
    if rc < 1e-14 || isnan(rc)
        lambda = 1e-10 * max(max(abs(diag(K_block))), 1e-20);
        K_block = K_block + lambda * eye(nsp);
    end

    K_block_inv = K_block \ eye(nsp);

    for ii = 1:nsp
        for jj = 1:nsp
            if abs(K_block_inv(ii, jj)) > 1e-30
                ptr = ptr + 1;
                I_idx(ptr) = idx(ii);
                J_idx(ptr) = idx(jj);
                V_val(ptr) = K_block_inv(ii, jj);
            end
        end
    end
end

% K_vv inverse (diagonal)
for si = 1:n_sensors
    ptr = ptr + 1;
    I_idx(ptr) = n_total + si;
    J_idx(ptr) = n_total + si;
    V_val(ptr) = 1.0 / max(K_vv_vec(si), 1e-30);
end

I_idx = I_idx(1:ptr);
J_idx = J_idx(1:ptr);
V_val = V_val(1:ptr);
n_aug = n_total + n_sensors;
W = sparse(I_idx, J_idx, V_val, n_aug, n_aug);
end

%% ========================================================================
%  DIAGONAL EXTRACTION: diag(inv(AtWA))
%  ========================================================================

function d = extract_diagonal_of_sparse_inverse_local(AtWA, batch_size)
% Extract diagonal of inv(AtWA) using sparse Cholesky + batch solves

n = size(AtWA, 1);
d = zeros(n, 1);

% Ensure positive definiteness
min_diag = min(diag(AtWA));
if min_diag <= 0
    lambda = abs(min_diag) + 1e-12 * max(abs(diag(AtWA)));
    AtWA = AtWA + lambda * speye(n);
end

% Use decomposition for repeated solves (caches factorization)
try
    dA = decomposition(AtWA, 'chol');
    use_decomp = true;
catch
    lambda = 1e-10 * max(abs(diag(AtWA)));
    AtWA = AtWA + lambda * speye(n);
    try
        dA = decomposition(AtWA, 'ldl');
        use_decomp = true;
    catch
        use_decomp = false;
    end
end

for start = 1:batch_size:n
    idx = start:min(start + batch_size - 1, n);
    nb = length(idx);

    E = sparse(idx, 1:nb, ones(nb, 1), n, nb);

    if use_decomp
        X = dA \ E;
    else
        X = AtWA \ E;
    end

    d(idx) = full(sum(E .* X, 1))';
end
end

%% ========================================================================
%  GLOBAL MEASUREMENT MATRIX C
%  ========================================================================

function [C, K_vv_vec] = build_global_C_matrix_local(sensor_config, cell_index_map, ...
    n_total, n_sensors, nsp)
% Build global sparse C matrix (n_sensors x n_total)

I_idx = zeros(n_sensors, 1);
J_idx = zeros(n_sensors, 1);
K_vv_vec = zeros(n_sensors, 1);

dir_to_state = [2, 3, 4, 5];
si = 0;

% Depth sensors → state 1 (h)
for di = 1:sensor_config.n_depth
    si = si + 1;
    r = sensor_config.depth_sensors(di, 1);
    c = sensor_config.depth_sensors(di, 2);
    ci = cell_index_map(r, c);
    if ci > 0
        I_idx(si) = si;
        J_idx(si) = (ci - 1) * nsp + 1;
        K_vv_vec(si) = sensor_config.depth_noise_var(di);
    end
end

% Flow sensors → states 2-5 (Q_l/Q_r/Q_u/Q_d)
for fi = 1:sensor_config.n_flow
    si = si + 1;
    r = sensor_config.flow_sensors(fi, 1);
    c = sensor_config.flow_sensors(fi, 2);
    dir = sensor_config.flow_sensors(fi, 3);
    ci = cell_index_map(r, c);
    if ci > 0
        I_idx(si) = si;
        state_offset = dir_to_state(dir);
        J_idx(si) = (ci - 1) * nsp + state_offset;
        K_vv_vec(si) = sensor_config.flow_noise_var(fi);
    end
end

% Infiltration sensors → state 6 (F)
for ii = 1:sensor_config.n_infil
    si = si + 1;
    r = sensor_config.infil_sensors(ii, 1);
    c = sensor_config.infil_sensors(ii, 2);
    ci = cell_index_map(r, c);
    if ci > 0
        I_idx(si) = si;
        J_idx(si) = (ci - 1) * nsp + 6;
        K_vv_vec(si) = sensor_config.infil_noise_var(ii);
    end
end

valid = I_idx > 0 & J_idx > 0;
C = sparse(I_idx(valid), J_idx(valid), ones(sum(valid), 1), n_sensors, n_total);
end

%% ========================================================================
%  PHYSICS HELPERS
%  ========================================================================

function [Q_l, Q_r, Q_u, Q_d] = cell_fluxes_global_local(h_grid, row, col, ...
    n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask)
% Manning CA fluxes for a cell

h = h_grid(row, col);
h_safe = max(h, 0);
h_reg = sqrt(h_safe^2 + h_smooth^2);
h_eff = max(h_reg - h_0, h_smooth);

wse = Z(row, col) + h_safe;
wse_left = wse; wse_right = wse; wse_up = wse; wse_down = wse;

[ny_g, nx_g] = size(Z);

if col > 1 && active_mask(row, col-1)
    wse_left = Z(row, col-1) + max(h_grid(row, col-1), 0);
end
if col < nx_g && active_mask(row, col+1)
    wse_right = Z(row, col+1) + max(h_grid(row, col+1), 0);
end
if row > 1 && active_mask(row-1, col)
    wse_up = Z(row-1, col) + max(h_grid(row-1, col), 0);
end
if row < ny_g && active_mask(row+1, col)
    wse_down = Z(row+1, col) + max(h_grid(row+1, col), 0);
end

S_left  = max((wse - wse_left) / dx, 0);
S_right = max((wse - wse_right) / dx, 0);
S_up    = max((wse - wse_up) / dy, 0);
S_down  = max((wse - wse_down) / dy, 0);

S_total = S_left + S_right + S_up + S_down;
S_max = max([S_left, S_right, S_up, S_down]);

v_max = (1 / n_cell) * h_eff^(2/3) * sqrt(max(S_max, 0));
v_max = min(v_max, v_limit);

flow_factor = min(1.0, h_safe / h_flow_threshold);
flow_factor = flow_factor * flow_factor * (3 - 2 * flow_factor);

if S_total > 1e-12
    w_l = S_left / S_total;  w_r = S_right / S_total;
    w_u = S_up / S_total;    w_d = S_down / S_total;
else
    w_l = 0; w_r = 0; w_u = 0; w_d = 0;
end

Q_l = flow_factor * v_max * h_eff * dy * w_l;
Q_r = flow_factor * v_max * h_eff * dy * w_r;
Q_u = flow_factor * v_max * h_eff * dx * w_u;
Q_d = flow_factor * v_max * h_eff * dx * w_d;
end

%% ========================================================================
%  UTILITY HELPERS
%  ========================================================================

function val = get_param_local(param, row, col)
if numel(param) > 1
    val = param(row, col);
else
    val = param;
end
end

function merged = merge_struct_defaults_local(base_struct, override_struct)
merged = base_struct;
fields = fieldnames(override_struct);
for ii = 1:numel(fields)
    merged.(fields{ii}) = override_struct.(fields{ii});
end
end

function R_val = get_rain_schedule_local(schedule, t)
R_val = 0;
for ri = 1:size(schedule, 1)
    if t >= schedule(ri, 1) && t < schedule(ri, 2)
        R_val = schedule(ri, 3);
        return;
    end
end
end

function [pattern_scalar, pattern_grid] = build_rain_pattern_local(base, ...
    use_compressed, n_active, compressed_to_grid, ny, nx, active_mask)
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

%% ========================================================================
%  TRAJECTORY EXTRACTION
%  ========================================================================

function [h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj] = ...
    extract_trajectories_global(x_out, nt, ny, nx, nsp, ...
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
            i = compressed_to_grid(c, 1);
            j = compressed_to_grid(c, 2);
            b = (c - 1) * nsp;
            h_traj(i,j,k) = x_k(b+1);
            Q_l_traj(i,j,k) = x_k(b+2);
            Q_r_traj(i,j,k) = x_k(b+3);
            Q_u_traj(i,j,k) = x_k(b+4);
            Q_d_traj(i,j,k) = x_k(b+5);
            F_traj(i,j,k) = x_k(b+6);
        end
    else
        for idx = 1:(nx*ny)
            i = ceil(idx / nx);
            j = mod(idx - 1, nx) + 1;
            b = (idx - 1) * nsp;
            h_traj(i,j,k) = x_k(b+1);
            Q_l_traj(i,j,k) = x_k(b+2);
            Q_r_traj(i,j,k) = x_k(b+3);
            Q_u_traj(i,j,k) = x_k(b+4);
            Q_d_traj(i,j,k) = x_k(b+5);
            F_traj(i,j,k) = x_k(b+6);
        end
    end
end
end

function y_mean = build_measurement_trajectory_global_local(h_traj, Q_l_traj, Q_r_traj, ...
    Q_u_traj, Q_d_traj, F_traj, sensor_config, cell_index_map, nt) %#ok<INUSL>
% Build mean measurement trajectory y = C * x

n_total = sensor_config.n_total;
y_mean = zeros(n_total, nt);
state_trajs = {h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, F_traj};

si = 0;
% Depth sensors
for di = 1:sensor_config.n_depth
    si = si + 1;
    r = sensor_config.depth_sensors(di, 1);
    c = sensor_config.depth_sensors(di, 2);
    y_mean(si, :) = squeeze(h_traj(r, c, :))';
end
% Flow sensors
dir_to_state = [2, 3, 4, 5];
for fi = 1:sensor_config.n_flow
    si = si + 1;
    r = sensor_config.flow_sensors(fi, 1);
    c = sensor_config.flow_sensors(fi, 2);
    dir = sensor_config.flow_sensors(fi, 3);
    traj = state_trajs{dir_to_state(dir)};
    y_mean(si, :) = squeeze(traj(r, c, :))';
end
% Infiltration sensors
for ii = 1:sensor_config.n_infil
    si = si + 1;
    r = sensor_config.infil_sensors(ii, 1);
    c = sensor_config.infil_sensors(ii, 2);
    y_mean(si, :) = squeeze(F_traj(r, c, :))';
end
end
