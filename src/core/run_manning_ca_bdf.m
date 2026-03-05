function [t_out, x_out, results] = run_manning_ca_bdf(params, R_input, duration, varargin)
% run_manning_ca_bdf: Solve Manning CA DAE using BDF time discretization
%
% TRUE DAE FORMULATION with 6 states per cell:
%   Differential: h (depth), F (cumulative infiltration)
%   Algebraic: Q_l, Q_r, Q_u, Q_d (discharges)
%
% DAE System:
%   dh/dt = R - f_infil + (Q_in - Q_out) / A       (differential)
%   dF/dt = f_infil                                 (differential)
%   0 = Q_l - Q_manning(h, S_l)                     (algebraic)
%   0 = Q_r - Q_manning(h, S_r)                     (algebraic)
%   0 = Q_u - Q_manning(h, S_u)                     (algebraic)
%   0 = Q_d - Q_manning(h, S_d)                     (algebraic)
%
% BDF1 (Backward Euler):
%   x_{n+1} - x_n = dt * f(x_{n+1})    for differential states
%   0 = g(x_{n+1})                      for algebraic states
%
% This is a proper Index-1 DAE solved with fixed-step BDF and Newton iteration.
% Unconditionally stable, works on steep terrain like Walnut Gulch.
%
% Inputs:
%   params    - Parameter structure (from setup_*_domain)
%   R_input   - Rainfall intensity [m/s] (scalar or vector per cell)
%   duration  - Simulation duration [seconds]
%
% Optional:
%   'dt'      - Time step (default: 1.0 s)
%   'BDFOrder' - 1 (Backward Euler) or 2 (BDF2) (default: 1)
%   'NewtonTol' - Newton iteration tolerance (default: 1e-8)
%   'MaxNewtonIter' - Max Newton iterations (default: 50)
%   'AdaptiveTimestep' - Enable adaptive dt reduction (default: true)
%   'MinTimeStep' - Minimum allowed timestep (default: 0.1 s)
%   'DisplayProgress' - Show progress (default: true)
%   'OutputInterval' - Save every N steps (default: 60)
%
% Outputs:
%   t_out   - Time vector [s]
%   x_out   - State matrix (n_times x n_states) - [h, Q_l, Q_r, Q_u, Q_d, F] per cell
%   results - Structure with water balance, computation time, etc.
%
% Author: M. H. Kazma
% Date: January 2026

%% Parse inputs
p = inputParser;
addParameter(p, 'dt', 1.0);
addParameter(p, 'BDFOrder', 1);
addParameter(p, 'NewtonTol', 1e-6);  % Reasonable tolerance for practical problems
addParameter(p, 'MaxNewtonIter', 50);
addParameter(p, 'AdaptiveTimestep', true);
addParameter(p, 'MinTimeStep', 0.1);
addParameter(p, 'DisplayProgress', true);
addParameter(p, 'OutputInterval', 60);
parse(p, varargin{:});

dt_base = p.Results.dt;
bdf_order = p.Results.BDFOrder;
newton_tol = p.Results.NewtonTol;
max_newton_iter = p.Results.MaxNewtonIter;
use_adaptive_dt = p.Results.AdaptiveTimestep;
min_dt = p.Results.MinTimeStep;
display_progress = p.Results.DisplayProgress;
output_interval = p.Results.OutputInterval;

%% Extract parameters
grid_size = params.grid_size;
nx = grid_size(2);
ny = grid_size(1);
dx = params.dx;
dy = params.dy;
cell_area = dx * dy;

Z = params.Z;
n_manning = params.n_manning;
h_0 = params.h_0;
F_min = params.F_min;
F_initial = F_min;  % Default: start at F_min
if isfield(params, 'F_initial')
    F_initial = max(params.F_initial, F_min);  % Must be >= F_min
end

% Compressed indexing
if isfield(params, 'use_compressed_indexing') && params.use_compressed_indexing
    use_compressed = true;
    n_cells = params.valid_cells;
    grid_to_compressed = params.grid_to_compressed;
    compressed_to_grid = params.compressed_to_grid;
else
    use_compressed = false;
    n_cells = nx * ny;
    grid_to_compressed = [];
    compressed_to_grid = [];
end

% 6 states per cell: [h, Q_l, Q_r, Q_u, Q_d, F]
states_per_cell = 6;
n_states = n_cells * states_per_cell;

% Prepare rainfall (can be scalar, vector, time-dependent struct, or rain_matrix)
use_rain_schedule = false;
use_rain_matrix = false;
if isstruct(R_input) && isfield(R_input, 'rain_matrix')
    % Spatially distributed, time-varying rainfall matrix [n_cells x n_times]
    rain_matrix = R_input.rain_matrix;
    dt_rain = R_input.dt_rain;  % Time resolution of matrix columns [seconds]
    use_rain_matrix = true;
    use_time_dependent_rain = true;
    R_base_vec = mean(rain_matrix, 2);  % Time-averaged per-cell rate for balance
    actual_rain_duration = size(rain_matrix, 2) * dt_rain;
    R_for_balance = R_base_vec;  % Approximate for water balance
    if display_progress
        fprintf('  Rain matrix: %d cells x %d timesteps (dt_rain=%.0fs)\n', ...
            size(rain_matrix, 1), size(rain_matrix, 2), dt_rain);
        fprintf('  Peak spatial intensity: %.1f mm/hr\n', ...
            max(rain_matrix(:)) * 1000 * 3600);
    end
elseif isstruct(R_input) && isfield(R_input, 'schedule')
    % Time-dependent rainfall schedule: [t_start, t_end, intensity]
    schedule = R_input.schedule;
    use_rain_schedule = true;
    use_time_dependent_rain = true;
    
    if isfield(R_input, 'pattern')
        R_pattern = R_input.pattern;
    else
        R_pattern = 1;
    end
    
    if isscalar(R_pattern)
        R_base_vec = R_pattern * ones(n_cells, 1);
    else
        R_base_vec = R_pattern(:);
    end
    
    % Store for water balance calculation
    actual_rain_duration = duration;
    R_for_balance = R_base_vec;
elseif isstruct(R_input)
    % Time-dependent rainfall with duration
    R_base = R_input.base;  % Base intensity (m/s)
    rain_duration = R_input.duration;  % Duration (seconds)
    use_time_dependent_rain = true;
    if isscalar(R_base)
        R_base_vec = R_base * ones(n_cells, 1);
    else
        R_base_vec = R_base(:);
    end
    % Store for water balance calculation
    actual_rain_duration = rain_duration;
    R_for_balance = R_base_vec;  % Original intensity for balance calc
elseif isscalar(R_input)
    R_base_vec = R_input * ones(n_cells, 1);
    use_time_dependent_rain = false;
    actual_rain_duration = duration;  % Full simulation duration
    R_for_balance = R_base_vec;
else
    R_base_vec = R_input(:);
    use_time_dependent_rain = false;
    actual_rain_duration = duration;  % Full simulation duration
    R_for_balance = R_base_vec;
end

% Infiltration parameters
if isfield(params, 'Ks')
    Ks = params.Ks;
    eta = params.eta;
    psi_f = params.psi_f;
    use_infiltration = true;
else
    use_infiltration = false;
    Ks = 0; eta = 0; psi_f = 0;
end

% Regularization parameters - read from params if available, else use defaults
if isfield(params, 'h_smooth')
    h_smooth = params.h_smooth;
else
    h_smooth = 5e-3;        % 5 mm smoothing for h^(2/3) (default for steep terrain)
end
if isfield(params, 'h_flow_threshold')
    h_flow_threshold = params.h_flow_threshold;
else
    h_flow_threshold = 2e-3; % 2 mm threshold for flow
end
v_limit = 8.0;          % Max velocity m/s

% Outlet boundary parameters
if isfield(params, 'outlet_direction')
    outlet_direction = params.outlet_direction;
else
    outlet_direction = 'down';
end
if isfield(params, 'slope_outlet')
    slope_outlet = params.slope_outlet;
else
    slope_outlet = 0.01;
end
if isfield(params, 'outlet_indices')
    outlet_indices = params.outlet_indices;
else
    outlet_indices = [];
end

% Pre-compute outlet mask for fast lookup (avoid ismember in loops)
is_outlet_cell = false(n_cells, 1);
is_outlet_cell(outlet_indices) = true;

%% Display setup
if display_progress
    fprintf('========================================\n');
    fprintf('Manning CA DAE - BDF Solver [v2024-01-27 FREE DRAINAGE FIX]\n');
    fprintf('========================================\n');
    fprintf('Formulation: TRUE DAE with 6 states per cell\n');
    fprintf('  Differential: h, F\n');
    fprintf('  Algebraic: Q_l, Q_r, Q_u, Q_d\n');
    fprintf('Domain: %d × %d = %d cells\n', ny, nx, n_cells);
    fprintf('States: %d (6 per cell)\n', n_states);
    fprintf('Duration: %.1f seconds (%.1f minutes)\n', duration, duration/60);
    fprintf('Time step: %.2f s (BDF%d)%s\n', dt_base, bdf_order, ...
        ternary(use_adaptive_dt, ' [adaptive]', ''));
    fprintf('Newton tolerance: %.0e\n', newton_tol);
    if use_adaptive_dt
        fprintf('Adaptive timestepping: enabled (min dt = %.2f s)\n', min_dt);
    end
    fprintf('\n*** USING TERRAIN-BASED FREE DRAINAGE AT BOUNDARIES ***\n\n');
end

%% Initialize state vector
n_steps = ceil(duration / dt_base);
n_output = ceil(n_steps / output_interval) + 1;

t_out = zeros(n_output, 1);
x_out = zeros(n_output, n_states);

% Initial state: dry surface, antecedent soil moisture
x = zeros(n_states, 1);
for c = 1:n_cells
    base = (c-1)*6;
    x(base + 1) = 0;          % h = 0 (dry surface)
    x(base + 2) = 0;          % Q_l = 0
    x(base + 3) = 0;          % Q_r = 0
    x(base + 4) = 0;          % Q_u = 0
    x(base + 5) = 0;          % Q_d = 0
    x(base + 6) = F_initial;  % F = F_initial (antecedent soil moisture)
end

% For BDF2
x_prev = x;

% Store initial
t_out(1) = 0;
x_out(1, :) = x';
Q_boundary_ts(1) = 0;  % No outflow initially
Q_outlet_ts = zeros(n_output, 1);  % Outlet-only (flume) discharge
Q_outlet_ts(1) = 0;
output_idx = 2;

if display_progress
    if F_initial > F_min
        fprintf('Initial: DRY surface, antecedent F=%.1f mm (F_min=%.1f mm)\n\n', F_initial*1000, F_min*1000);
    else
        fprintf('Initial: DRY (h=0, Q=0, F=F_min=%.1f mm)\n\n', F_min*1000);
    end
    fprintf('Running BDF DAE time-stepping...\n');
end

%% Time-stepping loop
tic;
dt = dt_base;  % Current timestep (may adapt)
n_timestep_reductions = 0;
n_timestep_increases = 0;
t = 0;
step = 0;
total_newton_iter = 0;
total_outflow = 0;
R_prev = zeros(n_cells, 1);  % AWR2026: Track previous R for rain-aware guess

while t < duration
    step = step + 1;
    t_new = min(t + dt, duration);
    dt_actual = t_new - t;
    
    % For Index-1 DAE with explicit algebraic constraints Q = g(h),
    % we solve using a semi-implicit approach:
    % 1. Compute Q from current h (explicit evaluation of algebraic constraint)
    % 2. Update h implicitly using BDF with fixed-point iteration
    
    % Apply time-dependent rainfall if configured
    if use_rain_matrix
        % Spatial IDW: look up column from pre-computed matrix
        t_idx = min(floor(t / dt_rain) + 1, size(rain_matrix, 2));
        R = rain_matrix(:, t_idx);
    elseif use_rain_schedule
        R_scalar = get_rainfall_from_schedule(schedule, t);
        R = R_scalar * R_base_vec;
    elseif use_time_dependent_rain
        if t < rain_duration
            R = R_base_vec;  % Rainfall active
        else
            if sum(R) > 0  % Only print once when transitioning
                fprintf('  *** RAINFALL STOPPED at t=%.1f s (%.1f min) ***\n', t, t/60);
            end
            R = zeros(n_cells, 1);  % Rainfall stopped
        end
    else
        R = R_base_vec;  % Constant rainfall throughout
    end
    
    x_new = x;  % Initial guess
    
    % AWR2026: Rain-aware predictor — adjust h for change in rainfall
    % When R drops (e.g., rain stops), h_{n+1} will be lower than h_n by
    % roughly dt*(R_old - R_new). Incorporating this into the initial guess
    % reduces Newton iterations at rain onset/cessation transitions.
    if step > 1
        dR = R - R_prev;  % Change in rainfall rate [m/s]
        for cc = 1:n_cells
            if dR(cc) ~= 0
                base_cc = (cc-1)*6;
                x_new(base_cc + 1) = max(x_new(base_cc + 1) + dt_actual * dR(cc), 0);
            end
        end
    end
    R_prev = R;  % Store for next step
    
    % First, update algebraic states based on current h
    x_new = update_algebraic_states(x_new, params, use_compressed, n_cells, ...
        grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, ...
        Z, n_manning, h_0, h_smooth, h_flow_threshold, v_limit, ...
        outlet_direction, slope_outlet, is_outlet_cell, dt_actual);
    
    % Under-relaxation for fixed-point convergence
    omega = 1.0;          % Relaxation factor (1.0 = full step)
    prev_res_norm = inf;
    
    %% AWR2026: CFL-based timestep ceiling for fixed-point stability
    if use_adaptive_dt
        h_max_cur = 0;
        for cc = 1:n_cells
            h_max_cur = max(h_max_cur, x((cc-1)*6 + 1));
        end
        if h_max_cur > 0.01  % Only when meaningful depth exists
            n_min_val = min(n_manning(n_manning(:) > 0));
            S_typical = 0.05;
            v_est = (1/n_min_val) * h_max_cur^(2/3) * sqrt(S_typical);
            dt_cfl = min(dx, dy) / max(v_est, 0.01);
            dt_cfl = max(dt_cfl, min_dt);
            if dt > dt_cfl
                dt = dt_cfl;
                t_new = min(t + dt, duration);
                dt_actual = t_new - t;
            end
        end
    end
    
    for newton_iter = 1:max_newton_iter
        % Compute the right-hand side dh/dt using current Q
        % This is a fixed-point iteration: h^{n+1} = h^n + dt * f(h^{n+1}, Q(h^{k}))
        
        [res_diff, Q_boundary, dh_dt_vec, dF_dt_vec] = compute_differential_residual_with_rhs(x_new, x, x_prev, ...
            dt_actual, step, bdf_order, R, params, ...
            use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
            nx, ny, dx, dy, cell_area, Z, n_manning, h_0, F_min, ...
            Ks, eta, psi_f, use_infiltration);
        
        res_norm = norm(res_diff);
        
        % Debug: print residual for first few steps if not converging
        if step <= 3 && newton_iter <= 10
            fprintf('      Newton iter %d: res_norm = %.6e (tol = %.6e)\n', newton_iter, res_norm, newton_tol);
            fprintf('        max(dh_dt) = %.6e, max(h_new) = %.6e\n', max(dh_dt_vec), max(x_new(1:6:end)));
            fprintf('        R(1) = %.6e, sample dh_dt(1) = %.6e\n', R(1), dh_dt_vec(1));
        end
        
        if res_norm < newton_tol
            break;
        end
        
        % AWR2026: Early exit for dt reduction — don't waste iterations
        if use_adaptive_dt && newton_iter > 15 && res_norm > newton_tol * 10 && dt > min_dt * 1.01
            break;  % Will trigger dt reduction below
        end
        
        % AWR2026: Adaptive under-relaxation from iteration 1 (was >2)
        % Floor lowered from 0.2 to 0.02 for h^(5/3) nonlinearity at large depths
        if newton_iter >= 1 && res_norm > prev_res_norm * 0.95
            omega = max(omega * 0.5, 0.02);  % Much lower floor for deep flows
        elseif newton_iter > 1 && res_norm < prev_res_norm * 0.5
            omega = min(omega * 1.2, 1.0);  % Recovering - relax damping
        end
        prev_res_norm = res_norm;
        
        % Fixed-point update with under-relaxation:
        % h_target = h_old + dt * dh_dt  (standard fixed-point)
        % h_new = (1-omega)*h_current + omega*h_target  (damped)
        
        for c = 1:n_cells
            base = (c-1)*6;
            
            if bdf_order == 1 || step == 1
                h_target = x(base + 1) + dt_actual * dh_dt_vec(c);
                F_target = x(base + 6) + dt_actual * dF_dt_vec(c);
            else
                h_target = (4/3)*x(base + 1) - (1/3)*x_prev(base + 1) + (2/3)*dt_actual * dh_dt_vec(c);
                F_target = (4/3)*x(base + 6) - (1/3)*x_prev(base + 6) + (2/3)*dt_actual * dF_dt_vec(c);
            end
            
            % Under-relaxed update
            x_new(base + 1) = (1 - omega) * x_new(base + 1) + omega * h_target;
            x_new(base + 6) = (1 - omega) * x_new(base + 6) + omega * F_target;
            
            % Enforce bounds
            x_new(base + 1) = max(x_new(base + 1), 0);
            x_new(base + 6) = max(x_new(base + 6), F_min);
        end
        
        % Update Q from new h (maintain algebraic constraint)
        x_new = update_algebraic_states(x_new, params, use_compressed, n_cells, ...
            grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, ...
            Z, n_manning, h_0, h_smooth, h_flow_threshold, v_limit, ...
            outlet_direction, slope_outlet, is_outlet_cell, dt_actual);
    end
    
    % AWR2026: Check for convergence failure and adapt timestep
    % Tightened acceptance threshold from 100x to 10x tolerance
    if use_adaptive_dt && (newton_iter >= max_newton_iter || ...
            (newton_iter > 15 && res_norm > newton_tol * 10))
        if res_norm > newton_tol * 10 && dt > min_dt * 1.01
            % Failed to converge — reduce timestep and retry
            dt = max(dt * 0.5, min_dt);
            n_timestep_reductions = n_timestep_reductions + 1;
            if display_progress && mod(n_timestep_reductions, 5) == 1
                fprintf('  WARNING: Newton not converging (res=%.2e, omega=%.2f, it=%d) at t=%.1f s. Reducing dt to %.4f s\n', ...
                    res_norm, omega, newton_iter, t, dt);
            end
            % Retry this timestep with smaller dt
            step = step - 1;
            continue;
        else
            % Residual is small enough — accept solution and continue
            if display_progress && mod(step, 60) == 0
                fprintf('  Note: Accepted with res=%.2e (omega=%.2f) at t=%.1f s\n', res_norm, omega, t);
            end
        end
    elseif use_adaptive_dt && newton_iter < 5 && dt < dt_base * 0.99
        % Converging easily - try increasing timestep
        dt = min(dt * 1.5, dt_base);
        n_timestep_increases = n_timestep_increases + 1;
    end
    
    % C
    total_newton_iter = total_newton_iter + newton_iter;
    
    % Compute boundary outflow from the FINAL converged state
    % This must be done AFTER the last update_algebraic_states call
    Q_boundary = compute_boundary_outflow_from_state(x_new, use_compressed, n_cells, ...
        compressed_to_grid, grid_to_compressed, nx, ny);
    
    % Debug: Print Q_boundary for first few timesteps
    if step <= 5
        fprintf('    DEBUG: Step %d, Q_boundary = %.6e m³/s\n', step, Q_boundary);
    end
    
    % Track boundary outflow
    total_outflow = total_outflow + Q_boundary * dt_actual;
    
    % Update states for next step
    x_prev = x;
    x = x_new;
    t = t_new;
    
    % Store output
    if mod(step, output_interval) == 0 || t >= duration
        t_out(output_idx) = t;
        x_out(output_idx, :) = x';
        Q_boundary_ts(output_idx) = Q_boundary;  % Save instantaneous boundary outflow

        % Compute Q_outlet: flows leaving domain at designated outlet cells only
        Q_outlet_now = compute_outlet_only_Q(x_new, outlet_indices, use_compressed, ...
            compressed_to_grid, grid_to_compressed, nx, ny);
        Q_outlet_ts(output_idx) = Q_outlet_now;

        output_idx = output_idx + 1;
        
        if display_progress
            h_max = 0;
            for c = 1:n_cells
                h_max = max(h_max, x((c-1)*6 + 1));
            end
            fprintf('  t = %6.1f s (%5.1f min), h_max = %.4f m, Newton = %d\n', ...
                t, t/60, h_max, newton_iter);
        end
    end
end

comp_time = toc;

% Trim output arrays
t_out = t_out(1:output_idx-1);
x_out = x_out(1:output_idx-1, :);
Q_boundary_ts = Q_boundary_ts(1:output_idx-1);
Q_outlet_ts = Q_outlet_ts(1:output_idx-1);

%% Extract final results
h_final = zeros(ny, nx);
F_final = zeros(ny, nx);
Q_out_final = zeros(ny, nx, 4);  % [Q_l, Q_r, Q_u, Q_d]

if use_compressed
    for c = 1:n_cells
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
        base = (c-1)*6;
        h_final(i, j) = x(base + 1);
        Q_out_final(i, j, 1) = x(base + 2);
        Q_out_final(i, j, 2) = x(base + 3);
        Q_out_final(i, j, 3) = x(base + 4);
        Q_out_final(i, j, 4) = x(base + 5);
        F_final(i, j) = x(base + 6);
    end
else
    for i = 1:ny
        for j = 1:nx
            c = (i-1)*nx + j;
            base = (c-1)*6;
            h_final(i, j) = x(base + 1);
            Q_out_final(i, j, 1) = x(base + 2);
            Q_out_final(i, j, 2) = x(base + 3);
            Q_out_final(i, j, 3) = x(base + 4);
            Q_out_final(i, j, 4) = x(base + 5);
            F_final(i, j) = x(base + 6);
        end
    end
end

%% Compute water balance
% Calculate total rainfall input (use actual rain duration, not simulation duration)
if use_rain_matrix
    % Exact: sum per-cell depth over all timesteps within simulation window
    n_rain_steps = min(size(rain_matrix, 2), floor(duration / dt_rain));
    if n_rain_steps > 0
        total_rainfall = sum(sum(rain_matrix(:, 1:n_rain_steps))) * dt_rain * cell_area;
    else
        total_rainfall = 0;
    end
    actual_rain_duration = n_rain_steps * dt_rain;
elseif use_rain_schedule
    total_rain_depth = 0;
    total_rain_time = 0;
    for rr = 1:size(schedule, 1)
        t0 = max(0, schedule(rr, 1));
        t1 = min(duration, schedule(rr, 2));
        if t1 > t0 && schedule(rr, 3) > 0
            total_rain_depth = total_rain_depth + schedule(rr, 3) * (t1 - t0);
            total_rain_time = total_rain_time + (t1 - t0);
        end
    end
    total_rainfall = total_rain_depth * sum(R_for_balance) * cell_area;
    actual_rain_duration = total_rain_time;
elseif length(R_for_balance) == n_cells
    % R is a vector: sum over all cells
    total_rainfall = sum(R_for_balance) * actual_rain_duration * cell_area;
else
    % R is scalar: multiply by number of cells
    total_rainfall = R_for_balance(1) * actual_rain_duration * n_cells * cell_area;
end

if use_compressed
    % For compressed indexing: only sum over valid cells
    valid_mask = params.valid_mask;
    total_ponding = sum(h_final(valid_mask)) * cell_area;
    
    if use_infiltration
        % Only count infiltration in valid cells (F - F_initial is the infiltrated depth)
        F_valid = F_final(valid_mask);
        total_infiltration = sum(max(F_valid - F_initial, 0)) * cell_area;
    else
        total_infiltration = 0;
    end
else
    % For full grid: sum all non-NaN cells
    h_valid = h_final(~isnan(h_final));
    total_ponding = sum(h_valid) * cell_area;
    
    if use_infiltration
        F_valid = F_final(~isnan(F_final));
        total_infiltration = sum(max(F_valid - F_initial, 0)) * cell_area;
    else
        total_infiltration = 0;
    end
end

% DEBUG: Check boundary outflow components
total_runoff = total_outflow;

% Additional mass balance diagnostics
total_water_initial = n_cells * F_initial * cell_area; % Initial F = F_initial everywhere
water_added = total_rainfall;
water_removed = total_runoff;

% Total water final = ponding + infiltrated water (not counting initial F_min)
total_water_final = total_ponding + total_infiltration;

% Expected change = input - output
expected_change = water_added - water_removed;
actual_change = total_water_final;
discrepancy = actual_change - expected_change;

% Count boundary cells and check their Q values
if use_compressed && display_progress
    n_boundary_cells = 0;
    total_boundary_Q = 0;
    Q_components = zeros(4, 1); % [left, right, up, down]
    
    for c = 1:n_cells
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
        base = (c-1)*6;
        Q_l = x(base + 2);
        Q_r = x(base + 3);
        Q_u = x(base + 4);
        Q_d = x(base + 5);
        
        % Check if this is a boundary cell
        is_boundary = false;
        if j == 1 || grid_to_compressed(i, j-1) == 0
            n_boundary_cells = n_boundary_cells + 1;
            total_boundary_Q = total_boundary_Q + Q_l;
            Q_components(1) = Q_components(1) + Q_l;
            is_boundary = true;
        end
        if j == nx || grid_to_compressed(i, j+1) == 0
            if ~is_boundary
                n_boundary_cells = n_boundary_cells + 1;
                is_boundary = true;
            end
            total_boundary_Q = total_boundary_Q + Q_r;
            Q_components(2) = Q_components(2) + Q_r;
        end
        if i == 1 || grid_to_compressed(i-1, j) == 0
            if ~is_boundary
                n_boundary_cells = n_boundary_cells + 1;
                is_boundary = true;
            end
            total_boundary_Q = total_boundary_Q + Q_u;
            Q_components(3) = Q_components(3) + Q_u;
        end
        if i == ny || grid_to_compressed(i+1, j) == 0
            if ~is_boundary
                n_boundary_cells = n_boundary_cells + 1;
            end
            total_boundary_Q = total_boundary_Q + Q_d;
            Q_components(4) = Q_components(4) + Q_d;
        end
    end
    if use_adaptive_dt
        fprintf('Timestep adaptations: %d reductions, %d increases\n', ...
            n_timestep_reductions, n_timestep_increases);
    end
end

balance_error = total_rainfall - total_ponding - total_infiltration - total_runoff;

%% Store results
results.h_final = h_final;
results.F_final = F_final;
results.Q_final = Q_out_final;
results.Q_boundary_ts = Q_boundary_ts;  % Total boundary outflow timeseries [m³/s]
results.Q_outlet_ts = Q_outlet_ts;      % Outlet (flume) cells only [m³/s]

% Store rainfall information for PSE (preserve struct format if time-dependent)
if use_rain_matrix
    results.R = struct('rain_matrix', rain_matrix, 'dt_rain', dt_rain);
elseif use_rain_schedule
    results.R = struct('schedule', schedule, 'pattern', R_base_vec);
elseif use_time_dependent_rain
    results.R = struct('base', mean(R_for_balance), 'duration', actual_rain_duration);
else
    results.R = mean(R_for_balance);  % Scalar for constant rainfall
end

results.rain_duration = actual_rain_duration;  % Store actual rain duration
results.computation_time = comp_time;
results.total_newton_iter = total_newton_iter;
results.avg_newton_per_step = total_newton_iter / step;

results.water_balance.rainfall = total_rainfall;
results.water_balance.ponding = total_ponding;
results.water_balance.infiltration = total_infiltration;
results.water_balance.runoff = total_runoff;
results.water_balance.error = balance_error;
results.water_balance.error_pct = 100 * abs(balance_error) / max(total_rainfall, 1);

if display_progress
    fprintf('\n========================================\n');
    fprintf('Simulation Complete\n');
    fprintf('========================================\n');
    fprintf('Computation time: %.2f s\n', comp_time);
    fprintf('Total Newton iterations: %d (avg %.1f/step)\n', ...
        total_newton_iter, results.avg_newton_per_step);
    fprintf('\nWater Balance:\n');
    fprintf('  Rainfall:      %10.1f m³ (100.0%%)\n', total_rainfall);
    fprintf('  Ponding:       %10.1f m³ (%5.1f%%)\n', total_ponding, 100*total_ponding/total_rainfall);
    fprintf('  Infiltration:  %10.1f m³ (%5.1f%%)\n', total_infiltration, 100*total_infiltration/total_rainfall);
    fprintf('  Runoff:        %10.1f m³ (%5.1f%%)\n', total_runoff, 100*total_runoff/total_rainfall);
    fprintf('  ----------------------------------------\n');
    fprintf('  Total Output:  %10.1f m³ (%5.1f%%)\n', total_ponding + total_infiltration + total_runoff, 100*(total_ponding + total_infiltration + total_runoff)/total_rainfall);
    fprintf('  Balance Error: %10.2f m³ (%5.3f%%)\n', balance_error, results.water_balance.error_pct);
    
    % Display boundary diagnostics
    if use_compressed && exist('n_boundary_cells', 'var')
        fprintf('\nBoundary Diagnostics:\n');
        fprintf('  Boundary cells: %d (%.1f%% of domain)\n', n_boundary_cells, 100*n_boundary_cells/n_cells);
        fprintf('  Final boundary Q (instantaneous): %.4e m³/s\n', total_boundary_Q);
        fprintf('    Left boundary:  %.4e m³/s\n', Q_components(1));
        fprintf('    Right boundary: %.4e m³/s\n', Q_components(2));
        fprintf('    Up boundary:    %.4e m³/s\n', Q_components(3));
        fprintf('    Down boundary:  %.4e m³/s\n', Q_components(4));
        fprintf('  Total boundary outflow (integrated): %.2f m³\n', total_outflow);
        
        fprintf('\nMass Balance Check:\n');
        fprintf('  Initial storage (F_initial): %10.1f m³\n', total_water_initial);
        fprintf('  Water input (rainfall):  %10.1f m³\n', water_added);
        fprintf('  Water output (runoff):   %10.1f m³\n', water_removed);
        fprintf('  Expected change:         %10.1f m³\n', expected_change);
        fprintf('  Actual change (h+F-F0):  %10.1f m³\n', actual_change);
        fprintf('  Discrepancy:             %10.2f m³ (%.4f%%)\n', discrepancy, 100*abs(discrepancy)/max(water_added, 1));
        
        fprintf('\n  Balance check: Rainfall = Ponding + Infiltration + Runoff + Error\n');
        fprintf('                 %.1f = %.1f + %.1f + %.1f + %.1f m³\n', ...
            total_rainfall, total_ponding, total_infiltration, total_runoff, balance_error);
    end
    
    % Display final state statistics (similar to semi-explicit method)
    fprintf('\nFinal State:\n');
    if use_compressed
        h_valid = h_final(valid_mask);
        F_valid = F_final(valid_mask);
        fprintf('  Water depth: min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
            min(h_valid)*1000, max(h_valid)*1000, mean(h_valid)*1000);
        fprintf('  Infiltration: min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
            min(F_valid)*1000, max(F_valid)*1000, mean(F_valid)*1000);
    else
        h_nonzero = h_final(h_final > 0);
        F_nonzero = F_final(F_final > F_min);
        if ~isempty(h_nonzero)
            fprintf('  Water depth: min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
                min(h_nonzero)*1000, max(h_nonzero)*1000, mean(h_nonzero)*1000);
        end
        if ~isempty(F_nonzero)
            fprintf('  Infiltration: min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
                min(F_nonzero)*1000, max(F_nonzero)*1000, mean(F_nonzero)*1000);
        end
    end
    
    %% Flow Direction Diagnostics
    fprintf('\nFlow Direction Diagnostics (final state):\n');
    Q_l_all = zeros(n_cells, 1); Q_r_all = zeros(n_cells, 1);
    Q_u_all = zeros(n_cells, 1); Q_d_all = zeros(n_cells, 1);
    Q_in_all = zeros(n_cells, 1);
    n_ponding = 0;
    infil_total_rate = 0;
    rain_total_rate = 0;
    for c = 1:n_cells
        base = (c-1)*6;
        Q_l_all(c) = x(base + 2);
        Q_r_all(c) = x(base + 3);
        Q_u_all(c) = x(base + 4);
        Q_d_all(c) = x(base + 5);
        if x(base + 1) > 1e-6
            n_ponding = n_ponding + 1;
        end
    end
    
    % Compute Q_in for each cell (same logic as residual)
    for c = 1:n_cells
        if use_compressed
            i = compressed_to_grid(c, 1);
            j = compressed_to_grid(c, 2);
        else
            i = floor((c-1) / nx) + 1;
            j = mod(c-1, nx) + 1;
        end
        Q_in_c = 0;
        if j > 1
            if use_compressed
                nc = grid_to_compressed(i, j-1);
            else
                nc = (i-1)*nx + (j-1);
            end
            if nc > 0, Q_in_c = Q_in_c + x((nc-1)*6 + 3); end  % Q_r from left
        end
        if j < nx
            if use_compressed
                nc = grid_to_compressed(i, j+1);
            else
                nc = (i-1)*nx + (j+1);
            end
            if nc > 0, Q_in_c = Q_in_c + x((nc-1)*6 + 2); end  % Q_l from right
        end
        if i > 1
            if use_compressed
                nc = grid_to_compressed(i-1, j);
            else
                nc = (i-2)*nx + j;
            end
            if nc > 0, Q_in_c = Q_in_c + x((nc-1)*6 + 5); end  % Q_d from upper
        end
        if i < ny
            if use_compressed
                nc = grid_to_compressed(i+1, j);
            else
                nc = i*nx + j;
            end
            if nc > 0, Q_in_c = Q_in_c + x((nc-1)*6 + 4); end  % Q_u from lower
        end
        Q_in_all(c) = Q_in_c;
    end
    
    Q_out_all = Q_l_all + Q_r_all + Q_u_all + Q_d_all;
    Q_net_all = Q_in_all - Q_out_all;
    
    fprintf('  Cells with ponding (h > 1e-6): %d / %d (%.1f%%)\n', n_ponding, n_cells, 100*n_ponding/n_cells);
    fprintf('  Q_l  (left):   nonzero=%d, sum=%.4f m³/s, max=%.4f m³/s\n', ...
        sum(Q_l_all > 1e-15), sum(Q_l_all), max(Q_l_all));
    fprintf('  Q_r  (right):  nonzero=%d, sum=%.4f m³/s, max=%.4f m³/s\n', ...
        sum(Q_r_all > 1e-15), sum(Q_r_all), max(Q_r_all));
    fprintf('  Q_u  (up/i-1): nonzero=%d, sum=%.4f m³/s, max=%.4f m³/s\n', ...
        sum(Q_u_all > 1e-15), sum(Q_u_all), max(Q_u_all));
    fprintf('  Q_d  (dn/i+1): nonzero=%d, sum=%.4f m³/s, max=%.4f m³/s\n', ...
        sum(Q_d_all > 1e-15), sum(Q_d_all), max(Q_d_all));
    fprintf('  Q_in  (total inflow):  sum=%.4f m³/s\n', sum(Q_in_all));
    fprintf('  Q_out (total outflow): sum=%.4f m³/s\n', sum(Q_out_all));
    fprintf('  Q_net (in-out):        sum=%.4f m³/s\n', sum(Q_net_all));
    fprintf('  Max single-cell Q_out: %.4f m³/s at cell %d\n', max(Q_out_all), find(Q_out_all == max(Q_out_all), 1));
    
    % Flow fraction by direction
    total_Q = sum(Q_l_all) + sum(Q_r_all) + sum(Q_u_all) + sum(Q_d_all);
    if total_Q > 0
        fprintf('  Flow fraction: L=%.1f%%, R=%.1f%%, U=%.1f%%, D=%.1f%%\n', ...
            100*sum(Q_l_all)/total_Q, 100*sum(Q_r_all)/total_Q, ...
            100*sum(Q_u_all)/total_Q, 100*sum(Q_d_all)/total_Q);
    end
    
    % Check mass balance at cell level
    % dh/dt should be ~ (R - f_infil + Q_net/A)
    % At end of sim with no rain, dh/dt should be negative (draining)
    fprintf('  Cells with Q_in > Q_out: %d (%.1f%%) — accumulating flow\n', ...
        sum(Q_net_all > 0), 100*sum(Q_net_all > 0)/n_cells);
    fprintf('  Cells with Q_in = Q_out = 0: %d (%.1f%%) — stagnant\n', ...
        sum(Q_in_all == 0 & Q_out_all == 0), 100*sum(Q_in_all == 0 & Q_out_all == 0)/n_cells);
    
    fprintf('========================================\n');
end

end

%% ========================================================================
%  COMPUTE DAE RESIDUAL
%  ========================================================================

function [res_diff, res_alg, Q_boundary] = compute_dae_residual(x, x_old, x_prev, ...
    dt, step, bdf_order, R, params, ...
    use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
    nx, ny, dx, dy, Z, n_manning, h_0, F_min, ...
    Ks, eta, psi_f, use_infiltration, ...
    h_smooth, h_flow_threshold, v_limit, ...
    outlet_direction, slope_outlet, outlet_indices)
% Compute residual for the DAE system
%
% Differential residual (for h, F):
%   BDF1: res = (x - x_old)/dt - f(x)
%   BDF2: res = (x - 4/3*x_old + 1/3*x_prev)/(2/3*dt) - f(x)
%
% Algebraic residual (for Q):
%   res = Q - Q_manning(h, S)

cell_area = dx * dy;

res_diff = zeros(n_cells * 2, 1);  % [h, F] residuals
res_alg = zeros(n_cells * 4, 1);   % [Q_l, Q_r, Q_u, Q_d] residuals
Q_boundary = 0;

for c = 1:n_cells
    % Get grid indices
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1) / nx) + 1;
        j = mod(c-1, nx) + 1;
    end
    
    base = (c-1)*6;
    
    % Extract current states
    h = max(x(base + 1), 0);
    Q_l = x(base + 2);
    Q_r = x(base + 3);
    Q_u = x(base + 4);
    Q_d = x(base + 5);
    F = max(x(base + 6), F_min);
    
    % Extract old states
    h_old = x_old(base + 1);
    F_old = x_old(base + 6);
    
    % For BDF2
    if bdf_order == 2 && step > 1
        h_prev = x_prev(base + 1);
        F_prev = x_prev(base + 6);
    end
    
    % Cell properties
    z = Z(i, j);
    n = n_manning(i, j);
    
    %% Compute infiltration rate
    if use_infiltration
        F_safe = max(F, F_min);
        if numel(Ks) > 1
            Ks_c = Ks(i, j);
            eta_c = eta(i, j);
            psi_c = psi_f(i, j);
        else
            Ks_c = Ks;
            eta_c = eta;
            psi_c = psi_f;
        end
        
        eta_h = eta_c * (h + psi_c);
        f_capacity = Ks_c * (1 + eta_h / F_safe);
        
        % Available water supply = ponded water from PREVIOUS timestep + rain
        % Using h_old (not current iteration h) prevents fixed-point oscillation.
        % When capacity > supply, all water infiltrates and h -> 0 in one iteration.
        % When capacity < supply, ponding accumulates at rate (R - f_capacity).
        i_a = max(h_old, 0) / dt + R(c);
        scale = max(f_capacity, i_a) + 1e-30;
        eps_smooth = 1e-12 * scale^2;
        f_infil = 0.5 * (f_capacity + i_a - sqrt((f_capacity - i_a)^2 + eps_smooth));
        f_infil = max(f_infil, 0);
    else
        f_infil = 0;
    end
    
    %% Compute inflows from neighbors
    Q_in = 0;
    
    % Left neighbor's Q_r flows into this cell
    if j > 1
        if use_compressed
            nc = grid_to_compressed(i, j-1);
        else
            nc = (i-1)*nx + (j-1);
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 3);  % Q_r from left neighbor
        end
    end
    
    % Right neighbor's Q_l flows into this cell
    if j < nx
        if use_compressed
            nc = grid_to_compressed(i, j+1);
        else
            nc = (i-1)*nx + (j+1);
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 2);  % Q_l from right neighbor
        end
    end
    
    % Upper neighbor's Q_d flows into this cell
    if i > 1
        if use_compressed
            nc = grid_to_compressed(i-1, j);
        else
            nc = (i-2)*nx + j;
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 5);  % Q_d from upper neighbor
        end
    end
    
    % Lower neighbor's Q_u flows into this cell
    if i < ny
        if use_compressed
            nc = grid_to_compressed(i+1, j);
        else
            nc = i*nx + j;
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 4);  % Q_u from lower neighbor
        end
    end
    
    Q_out = Q_l + Q_r + Q_u + Q_d;
    
    %% Differential residuals
    % dh/dt = R - f_infil + (Q_in - Q_out) / A
    dh_dt = R(c) - f_infil + (Q_in - Q_out) / cell_area;
    
    % dF/dt = f_infil
    dF_dt = f_infil;
    
    if bdf_order == 1 || step == 1
        % BDF1: (x - x_old)/dt = f(x)
        res_diff((c-1)*2 + 1) = (h - h_old) / dt - dh_dt;
        res_diff((c-1)*2 + 2) = (F - F_old) / dt - dF_dt;
    else
        % BDF2: (x - 4/3*x_old + 1/3*x_prev) / (2/3*dt) = f(x)
        res_diff((c-1)*2 + 1) = (h - (4/3)*h_old + (1/3)*h_prev) / ((2/3)*dt) - dh_dt;
        res_diff((c-1)*2 + 2) = (F - (4/3)*F_old + (1/3)*F_prev) / ((2/3)*dt) - dF_dt;
    end
    
    %% Algebraic residuals: 0 = Q - Q_manning(h, S)
    % Compute Manning Q for each direction
    
    h_safe = max(h, 0);
    h_reg = sqrt(h_safe^2 + h_smooth^2);
    h_eff = max(h_reg - h_0, h_smooth);
    
    wse = z + h_safe;
    
    % Get neighbor water surface elevations
    wse_l = wse; wse_r = wse; wse_u = wse; wse_d = wse;
    
    % Check if this is an outlet cell
    is_outlet = ismember(c, outlet_indices);
    
    % Left neighbor
    if j > 1
        if use_compressed
            nc = grid_to_compressed(i, j-1);
        else
            nc = (i-1)*nx + (j-1);
        end
        if nc > 0
            h_n = max(x((nc-1)*6 + 1), 0);
            wse_l = Z(i, j-1) + h_n;
        else
            % No valid neighbor: use terrain slope for free drainage
            terrain_slope = (z - Z(i, j-1)) / dx;
            if is_outlet
                wse_l = wse - slope_outlet * dx;
            else
                % Extrapolate using terrain slope (min slope for drainage)
                wse_l = wse - max(terrain_slope, 0.001) * dx;
            end
        end
    else
        % Boundary edge
        if strcmp(outlet_direction, 'left') || is_outlet
            wse_l = wse - slope_outlet * dx;
        else
            % Free drainage with minimum slope
            wse_l = wse - 0.001 * dx;
        end
    end
    
    % Right neighbor
    if j < nx
        if use_compressed
            nc = grid_to_compressed(i, j+1);
        else
            nc = (i-1)*nx + (j+1);
        end
        if nc > 0
            h_n = max(x((nc-1)*6 + 1), 0);
            wse_r = Z(i, j+1) + h_n;
        else
            % No valid neighbor: use terrain slope for free drainage
            terrain_slope = (z - Z(i, j+1)) / dx;
            if is_outlet
                wse_r = wse - slope_outlet * dx;
            else
                % Extrapolate using terrain slope (min slope for drainage)
                wse_r = wse - max(terrain_slope, 0.001) * dx;
            end
        end
    else
        % Boundary edge
        if strcmp(outlet_direction, 'right') || is_outlet
            wse_r = wse - slope_outlet * dx;
        else
            % Free drainage with minimum slope
            wse_r = wse - 0.001 * dx;
        end
    end
    
    % Upper neighbor
    if i > 1
        if use_compressed
            nc = grid_to_compressed(i-1, j);
        else
            nc = (i-2)*nx + j;
        end
        if nc > 0
            h_n = max(x((nc-1)*6 + 1), 0);
            wse_u = Z(i-1, j) + h_n;
        else
            % No valid neighbor: use terrain slope for free drainage
            terrain_slope = (z - Z(i-1, j)) / dy;
            if is_outlet
                wse_u = wse - slope_outlet * dy;
            else
                % Extrapolate using terrain slope (min slope for drainage)
                wse_u = wse - max(terrain_slope, 0.001) * dy;
            end
        end
    else
        % Boundary edge
        if strcmp(outlet_direction, 'up') || is_outlet
            wse_u = wse - slope_outlet * dy;
        else
            % Free drainage with minimum slope
            wse_u = wse - 0.001 * dy;
        end
    end
    
    % Lower neighbor
    if i < ny
        if use_compressed
            nc = grid_to_compressed(i+1, j);
        else
            nc = i*nx + j;
        end
        if nc > 0
            h_n = max(x((nc-1)*6 + 1), 0);
            wse_d = Z(i+1, j) + h_n;
        else
            % No valid neighbor: use terrain slope for free drainage
            terrain_slope = (z - Z(i+1, j)) / dy;
            if is_outlet
                wse_d = wse - slope_outlet * dy;
            else
                % Extrapolate using terrain slope (min slope for drainage)
                wse_d = wse - max(terrain_slope, 0.001) * dy;
            end
        end
    else
        % Boundary edge
        if strcmp(outlet_direction, 'down') || is_outlet
            wse_d = wse - slope_outlet * dy;
        else
            % Free drainage with minimum slope
            wse_d = wse - 0.001 * dy;
        end
    end
    
    % Compute slopes
    S_l = max((wse - wse_l) / dx, 0);
    S_r = max((wse - wse_r) / dx, 0);
    S_u = max((wse - wse_u) / dy, 0);
    S_d = max((wse - wse_d) / dy, 0);
    
    S_max = max([S_l, S_r, S_u, S_d]);
    S_total = S_l + S_r + S_u + S_d;
    
    % Velocity from Manning
    v_max = min((1/n) * h_eff^(2/3) * sqrt(max(S_max, 0)), v_limit);
    
    % Flow factor (smooth Hermite)
    ff = min(1.0, h_safe / h_flow_threshold);
    ff = ff * ff * (3 - 2*ff);
    
    % Compute Q_manning for each direction
    if S_total > 1e-12
        Q_manning_l = ff * v_max * h_eff * dy * (S_l / S_total);
        Q_manning_r = ff * v_max * h_eff * dy * (S_r / S_total);
        Q_manning_u = ff * v_max * h_eff * dx * (S_u / S_total);
        Q_manning_d = ff * v_max * h_eff * dx * (S_d / S_total);
    else
        Q_manning_l = 0;
        Q_manning_r = 0;
        Q_manning_u = 0;
        Q_manning_d = 0;
    end
    
    % Algebraic residuals: 0 = Q - Q_manning
    res_alg((c-1)*4 + 1) = Q_l - Q_manning_l;
    res_alg((c-1)*4 + 2) = Q_r - Q_manning_r;
    res_alg((c-1)*4 + 3) = Q_u - Q_manning_u;
    res_alg((c-1)*4 + 4) = Q_d - Q_manning_d;
    
    %% Track boundary outflow
    % Flow leaving the domain
    if j == 1 || (use_compressed && j > 1 && grid_to_compressed(i, j-1) == 0)
        Q_boundary = Q_boundary + Q_l;
    end
    if j == nx || (use_compressed && j < nx && grid_to_compressed(i, j+1) == 0)
        Q_boundary = Q_boundary + Q_r;
    end
    if i == 1 || (use_compressed && i > 1 && grid_to_compressed(i-1, j) == 0)
        Q_boundary = Q_boundary + Q_u;
    end
    if i == ny || (use_compressed && i < ny && grid_to_compressed(i+1, j) == 0)
        Q_boundary = Q_boundary + Q_d;
    end
end

end

%% ========================================================================
%  UPDATE ALGEBRAIC STATES (Direct solve for Q)
%  ========================================================================

function x = update_algebraic_states(x, params, use_compressed, n_cells, ...
    grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, ...
    Z, n_manning, h_0, h_smooth, h_flow_threshold, v_limit, ...
    outlet_direction, slope_outlet, is_outlet_cell, dt_step)
% Directly solve algebraic constraints: Q = Q_manning(h, S)
% This is valid because the algebraic constraint is explicit in Q
% Returns modified x with Q states updated
% is_outlet_cell is a pre-computed boolean array for fast lookup
% dt_step: current timestep [s] — used for flux limiter to prevent early-time overshoot

cell_area = dx * dy;  % [m^2]

for c = 1:n_cells
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1) / nx) + 1;
        j = mod(c-1, nx) + 1;
    end
    
    base = (c-1)*6;
    h = max(x(base + 1), 0);
    
    z = Z(i, j);
    n = n_manning(i, j);
    
    h_safe = max(h, 0);
    h_reg = sqrt(h_safe^2 + h_smooth^2);
    h_eff = max(h_reg - h_0, h_smooth);
    
    wse = z + h_safe;
    
    % Get neighbor WSEs (same logic as residual computation)
    % For boundary cells (no valid neighbor), allow free drainage using local slope
    wse_l = wse; wse_r = wse; wse_u = wse; wse_d = wse;
    is_outlet = is_outlet_cell(c);  % Fast boolean lookup instead of ismember
    
    % Default boundary slope for free drainage at watershed edges
    default_boundary_slope = 0.01;  % 1% slope for free drainage
    
    if j > 1
        if use_compressed
            nc = grid_to_compressed(i, j-1);
            if nc > 0
                wse_l = Z(i, j-1) + max(x((nc-1)*6 + 1), 0);
            else
                % No valid neighbor - apply free drainage BC
                % Allow water to leave watershed at boundary cells
                if is_outlet
                    wse_l = wse - slope_outlet * dx;
                elseif isnan(Z(i, j-1))
                    % Neighbor outside watershed - use default drainage slope
                    wse_l = wse - default_boundary_slope * dx;
                else
                    % Neighbor inside grid but masked - use terrain slope if downhill
                    local_slope = (z - Z(i, j-1)) / dx;
                    if local_slope > 0
                        wse_l = wse - local_slope * dx;
                    end
                end
            end
        else
            nc = (i-1)*nx + (j-1);
            wse_l = Z(i, j-1) + max(x((nc-1)*6 + 1), 0);
        end
    elseif strcmp(outlet_direction, 'left') || is_outlet
        wse_l = wse - slope_outlet * dx;
    end
    
    if j < nx
        if use_compressed
            nc = grid_to_compressed(i, j+1);
            if nc > 0
                wse_r = Z(i, j+1) + max(x((nc-1)*6 + 1), 0);
            else
                % No valid neighbor - apply free drainage BC
                if is_outlet
                    wse_r = wse - slope_outlet * dx;
                elseif isnan(Z(i, j+1))
                    wse_r = wse - default_boundary_slope * dx;
                else
                    local_slope = (z - Z(i, j+1)) / dx;
                    if local_slope > 0
                        wse_r = wse - local_slope * dx;
                    end
                end
            end
        else
            nc = (i-1)*nx + (j+1);
            wse_r = Z(i, j+1) + max(x((nc-1)*6 + 1), 0);
        end
    elseif strcmp(outlet_direction, 'right') || is_outlet
        wse_r = wse - slope_outlet * dx;
    end
    
    if i > 1
        if use_compressed
            nc = grid_to_compressed(i-1, j);
            if nc > 0
                wse_u = Z(i-1, j) + max(x((nc-1)*6 + 1), 0);
            else
                % No valid neighbor - apply free drainage BC
                if is_outlet
                    wse_u = wse - slope_outlet * dy;
                elseif isnan(Z(i-1, j))
                    wse_u = wse - default_boundary_slope * dy;
                else
                    local_slope = (z - Z(i-1, j)) / dy;
                    if local_slope > 0
                        wse_u = wse - local_slope * dy;
                    end
                end
            end
        else
            nc = (i-2)*nx + j;
            wse_u = Z(i-1, j) + max(x((nc-1)*6 + 1), 0);
        end
    elseif strcmp(outlet_direction, 'up') || is_outlet
        wse_u = wse - slope_outlet * dy;
    end
    
    if i < ny
        if use_compressed
            nc = grid_to_compressed(i+1, j);
            if nc > 0
                wse_d = Z(i+1, j) + max(x((nc-1)*6 + 1), 0);
            else
                % No valid neighbor - apply free drainage BC
                if is_outlet
                    wse_d = wse - slope_outlet * dy;
                elseif isnan(Z(i+1, j))
                    wse_d = wse - default_boundary_slope * dy;
                else
                    local_slope = (z - Z(i+1, j)) / dy;
                    if local_slope > 0
                        wse_d = wse - local_slope * dy;
                    end
                end
            end
        else
            nc = i*nx + j;
            wse_d = Z(i+1, j) + max(x((nc-1)*6 + 1), 0);
        end
    elseif strcmp(outlet_direction, 'down') || is_outlet
        % Bottom boundary - apply outlet BC
        wse_d = wse - slope_outlet * dy;
    end
    
    S_l = max((wse - wse_l) / dx, 0);
    S_r = max((wse - wse_r) / dx, 0);
    S_u = max((wse - wse_u) / dy, 0);
    S_d = max((wse - wse_d) / dy, 0);
    
    S_max = max([S_l, S_r, S_u, S_d]);
    S_total = S_l + S_r + S_u + S_d;
    
    v_max = min((1/n) * h_eff^(2/3) * sqrt(max(S_max, 0)), v_limit);
    
    ff = min(1.0, h_safe / h_flow_threshold);
    ff = ff * ff * (3 - 2*ff);
    
    if S_total > 1e-12
        x(base + 2) = ff * v_max * h_eff * dy * (S_l / S_total);
        x(base + 3) = ff * v_max * h_eff * dy * (S_r / S_total);
        x(base + 4) = ff * v_max * h_eff * dx * (S_u / S_total);
        x(base + 5) = ff * v_max * h_eff * dx * (S_d / S_total);
    else
        x(base + 2) = 0;
        x(base + 3) = 0;
        x(base + 4) = 0;
        x(base + 5) = 0;
    end
    
    %% AWR2026: Flux limiter — prevent draining more water than available
    % Matches Semi-Implicit approach (line 396): budget = h/dt + R + q_in
    % Hard linear clamp, same as CA_Routing and Semi-Implicit.
    if dt_step > 0
        Q_out_total = max(x(base+2),0) + max(x(base+3),0) + max(x(base+4),0) + max(x(base+5),0);
        % Incoming Q from neighbors (negative Q means inflow to this cell)
        Q_in_total = max(-x(base+2),0) + max(-x(base+3),0) + max(-x(base+4),0) + max(-x(base+5),0);
        available_vol = h * cell_area + Q_in_total * dt_step;  % [m^3]
        if Q_out_total * dt_step > available_vol && Q_out_total > 0
            scale = available_vol / (Q_out_total * dt_step);
            scale = max(scale, 0);  % safety
            for qs = 2:5
                if x(base + qs) > 0
                    x(base + qs) = x(base + qs) * scale;
                end
            end
        end
    end
end

end

%% ========================================================================
%  COMPUTE DIFFERENTIAL RESIDUAL (Simplified - Q substituted from constraint)
%  ========================================================================

function [res_diff, Q_boundary] = compute_differential_residual(x, x_old, x_prev, ...
    dt, step, bdf_order, R, params, ...
    use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
    nx, ny, dx, dy, cell_area, Z, n_manning, h_0, F_min, ...
    Ks, eta, psi_f, use_infiltration)
% Compute residual for differential equations only
% Q values are taken from current state (already updated via algebraic constraint)
%
% BDF1: res = (x - x_old)/dt - f(x)
% BDF2: res = (x - 4/3*x_old + 1/3*x_prev)/(2/3*dt) - f(x)
%
% USES SMOOTH MIN for infiltration to prevent oscillation

res_diff = zeros(n_cells * 2, 1);
Q_boundary = 0;

for c = 1:n_cells
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1) / nx) + 1;
        j = mod(c-1, nx) + 1;
    end
    
    base = (c-1)*6;
    
    % Current states
    h = max(x(base + 1), 0);
    Q_l = x(base + 2);
    Q_r = x(base + 3);
    Q_u = x(base + 4);
    Q_d = x(base + 5);
    F = max(x(base + 6), F_min);
    
    % Old states
    h_old = x_old(base + 1);
    F_old = x_old(base + 6);
    
    % For BDF2
    if bdf_order == 2 && step > 1
        h_prev = x_prev(base + 1);
        F_prev = x_prev(base + 6);
    end
    
    %% Infiltration
    if use_infiltration
        F_safe = max(F, F_min);
        if numel(Ks) > 1
            Ks_c = Ks(i, j);
            eta_c = eta(i, j);
            psi_c = psi_f(i, j);
        else
            Ks_c = Ks;
            eta_c = eta;
            psi_c = psi_f;
        end
        
        eta_h = eta_c * (h + psi_c);
        f_capacity = Ks_c * (1 + eta_h / F_safe);
        
        % Available water: ponded from previous timestep + current rain
        % Using h_old (not current iteration h) to avoid fixed-point oscillation
        i_a = max(h_old, 0) / dt + R(c);
        scale = max(f_capacity, i_a) + 1e-30;
        eps_smooth = 1e-12 * scale^2;
        f_infil = 0.5 * (f_capacity + i_a - sqrt((f_capacity - i_a)^2 + eps_smooth));
        f_infil = max(f_infil, 0);
    else
        f_infil = 0;
    end
    
    %% Inflows from neighbors
    Q_in = 0;
    
    if j > 1
        if use_compressed
            nc = grid_to_compressed(i, j-1);
        else
            nc = (i-1)*nx + (j-1);
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 3);  % Q_r from left
        end
    end
    
    if j < nx
        if use_compressed
            nc = grid_to_compressed(i, j+1);
        else
            nc = (i-1)*nx + (j+1);
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 2);  % Q_l from right
        end
    end
    
    if i > 1
        if use_compressed
            nc = grid_to_compressed(i-1, j);
        else
            nc = (i-2)*nx + j;
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 5);  % Q_d from upper
        end
    end
    
    if i < ny
        if use_compressed
            nc = grid_to_compressed(i+1, j);
        else
            nc = i*nx + j;
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 4);  % Q_u from lower
        end
    end
    
    Q_out = Q_l + Q_r + Q_u + Q_d;
    
    %% Right-hand side
    dh_dt = R(c) - f_infil + (Q_in - Q_out) / cell_area;
    dF_dt = f_infil;
    
    %% Residual
    if bdf_order == 1 || step == 1
        res_diff((c-1)*2 + 1) = (h - h_old) / dt - dh_dt;
        res_diff((c-1)*2 + 2) = (F - F_old) / dt - dF_dt;
    else
        res_diff((c-1)*2 + 1) = (h - (4/3)*h_old + (1/3)*h_prev) / ((2/3)*dt) - dh_dt;
        res_diff((c-1)*2 + 2) = (F - (4/3)*F_old + (1/3)*F_prev) / ((2/3)*dt) - dF_dt;
    end
end

end

%% ========================================================================
%  COMPUTE DIFFERENTIAL RESIDUAL WITH RHS VECTORS
%  ========================================================================

function [res_diff, Q_boundary, dh_dt_vec, dF_dt_vec] = compute_differential_residual_with_rhs(x, x_old, x_prev, ...
    dt, step, bdf_order, R, params, ...
    use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
    nx, ny, dx, dy, cell_area, Z, n_manning, h_0, F_min, ...
    Ks, eta, psi_f, use_infiltration)
% Compute residual for differential equations and also return dh/dt, dF/dt vectors
% for use in fixed-point iteration
% Also computes boundary outflow CONSISTENTLY with the mass balance
%
% USES SMOOTH MIN for infiltration to prevent oscillation

res_diff = zeros(n_cells * 2, 1);
dh_dt_vec = zeros(n_cells, 1);
dF_dt_vec = zeros(n_cells, 1);
Q_boundary = 0;

for c = 1:n_cells
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1) / nx) + 1;
        j = mod(c-1, nx) + 1;
    end
    
    base = (c-1)*6;
    
    % Current states
    h = max(x(base + 1), 0);
    Q_l = x(base + 2);
    Q_r = x(base + 3);
    Q_u = x(base + 4);
    Q_d = x(base + 5);
    F = max(x(base + 6), F_min);
    
    % Old states
    h_old = x_old(base + 1);
    F_old = x_old(base + 6);
    
    % For BDF2
    if bdf_order == 2 && step > 1
        h_prev = x_prev(base + 1);
        F_prev = x_prev(base + 6);
    end
    
    %% AWR2026: Compute Q_out and Q_in BEFORE infiltration
    %  so net flux can be included in the infiltration availability budget,
    %  matching the semi-implicit approach (h_available = h + dt*R + dh_flux).
    Q_out = Q_l + Q_r + Q_u + Q_d;
    
    Q_in = 0;
    if j > 1
        if use_compressed
            nc = grid_to_compressed(i, j-1);
        else
            nc = (i-1)*nx + (j-1);
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 3);  % Q_r from left
        end
    end
    if j < nx
        if use_compressed
            nc = grid_to_compressed(i, j+1);
        else
            nc = (i-1)*nx + (j+1);
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 2);  % Q_l from right
        end
    end
    if i > 1
        if use_compressed
            nc = grid_to_compressed(i-1, j);
        else
            nc = (i-2)*nx + j;
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 5);  % Q_d from upper
        end
    end
    if i < ny
        if use_compressed
            nc = grid_to_compressed(i+1, j);
        else
            nc = i*nx + j;
        end
        if nc > 0
            Q_in = Q_in + x((nc-1)*6 + 4);  % Q_u from lower
        end
    end
    
    %% Infiltration (with flux-aware availability)
    if use_infiltration
        F_safe = max(F, F_min);
        if numel(Ks) > 1
            Ks_c = Ks(i, j);
            eta_c = eta(i, j);
            psi_c = psi_f(i, j);
        else
            Ks_c = Ks;
            eta_c = eta;
            psi_c = psi_f;
        end
        
        eta_h = eta_c * (h + psi_c);
        f_capacity = Ks_c * (1 + eta_h / F_safe);
        
        % AWR2026: Available water for infiltration.
        % The DAE routes water out faster than explicit methods (Q instantly
        % reaches Manning steady-state), leaving h_old depleted and reducing
        % available water for infiltration — causing a ~6% infiltration deficit.
        %
        % Fix: use max(h, h_old) instead of h_old alone. During the rising
        % limb, h (current iterate) > h_old, giving infiltration access to
        % the water that's actually present. During recession, max = h_old
        % (no change). This can't cause downward oscillation because the
        % floor is h_old — i_a can only increase, never decrease, vs. the
        % original formula.
        h_for_avail = max(h, max(h_old, 0));
        Q_net_rate = (Q_in - Q_out) / cell_area;  % [m/s]
        i_a = h_for_avail / dt + R(c) + max(Q_net_rate, 0);
        scale = max(f_capacity, i_a) + 1e-30;
        eps_smooth = 1e-12 * scale^2;
        f_infil = 0.5 * (f_capacity + i_a - sqrt((f_capacity - i_a)^2 + eps_smooth));
        f_infil = max(f_infil, 0);
    else
        f_infil = 0;
    end
    
    %% Compute boundary outflow (flows leaving to invalid neighbors)
    % This MUST match the flows used in dh_dt for mass conservation
    if use_compressed
        % Left boundary: no valid neighbor to the left
        if j == 1 || grid_to_compressed(i, j-1) == 0
            Q_boundary = Q_boundary + Q_l;
        end
        % Right boundary
        if j == nx || grid_to_compressed(i, j+1) == 0
            Q_boundary = Q_boundary + Q_r;
        end
        % Top boundary
        if i == 1 || grid_to_compressed(i-1, j) == 0
            Q_boundary = Q_boundary + Q_u;
        end
        % Bottom boundary
        if i == ny || grid_to_compressed(i+1, j) == 0
            Q_boundary = Q_boundary + Q_d;
        end
    else
        if j == 1, Q_boundary = Q_boundary + Q_l; end
        if j == nx, Q_boundary = Q_boundary + Q_r; end
        if i == 1, Q_boundary = Q_boundary + Q_u; end
        if i == ny, Q_boundary = Q_boundary + Q_d; end
    end
    
    %% Right-hand side
    dh_dt = R(c) - f_infil + (Q_in - Q_out) / cell_area;
    dF_dt = f_infil;
    
    % Store for fixed-point iteration
    dh_dt_vec(c) = dh_dt;
    dF_dt_vec(c) = dF_dt;
    
    %% Residual
    if bdf_order == 1 || step == 1
        res_diff((c-1)*2 + 1) = (h - h_old) / dt - dh_dt;
        res_diff((c-1)*2 + 2) = (F - F_old) / dt - dF_dt;
    else
        res_diff((c-1)*2 + 1) = (h - (4/3)*h_old + (1/3)*h_prev) / ((2/3)*dt) - dh_dt;
        res_diff((c-1)*2 + 2) = (F - (4/3)*F_old + (1/3)*F_prev) / ((2/3)*dt) - dF_dt;
    end
end

end

%% ========================================================================
%  COMPUTE OUTLET-ONLY Q (flume cells)
%  ========================================================================

function Q_outlet = compute_outlet_only_Q(x, outlet_indices, use_compressed, ...
    compressed_to_grid, grid_to_compressed, nx, ny)
% Sum Q flowing out of domain ONLY at designated outlet (flume) cells.
% For vtilted (outlet_direction='down'), outlet_indices = entire bottom row,
% so Q_outlet == Q_boundary. For Walnut Gulch, outlet_indices = 4 flume cells.

Q_outlet = 0;
for k = 1:length(outlet_indices)
    c = outlet_indices(k);
    base = (c-1)*6;
    
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
        if j == 1 || grid_to_compressed(i, j-1) == 0
            Q_outlet = Q_outlet + x(base + 2);  % Q_l
        end
        if j == nx || grid_to_compressed(i, j+1) == 0
            Q_outlet = Q_outlet + x(base + 3);  % Q_r
        end
        if i == 1 || grid_to_compressed(i-1, j) == 0
            Q_outlet = Q_outlet + x(base + 4);  % Q_u
        end
        if i == ny || grid_to_compressed(i+1, j) == 0
            Q_outlet = Q_outlet + x(base + 5);  % Q_d
        end
    else
        i = floor((c-1) / nx) + 1;
        j = mod(c-1, nx) + 1;
        if j == 1,  Q_outlet = Q_outlet + x(base + 2); end
        if j == nx, Q_outlet = Q_outlet + x(base + 3); end
        if i == 1,  Q_outlet = Q_outlet + x(base + 4); end
        if i == ny, Q_outlet = Q_outlet + x(base + 5); end
    end
end
end

%% ========================================================================
%  COMPUTE BOUNDARY OUTFLOW
%  ========================================================================

function Q_boundary = compute_boundary_outflow(x, use_compressed, n_cells, ...
    grid_to_compressed, nx, ny, params)
% Sum Q flowing out of domain boundaries
% For compressed indexing, check each direction for missing neighbors (boundary)

Q_boundary = 0;

if use_compressed
    % For compressed indexing, we need compressed_to_grid to get (i,j)
    compressed_to_grid = params.compressed_to_grid;
    
    for c = 1:n_cells
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
        
        base = (c-1)*6;
        Q_l = x(base + 2);
        Q_r = x(base + 3);
        Q_u = x(base + 4);
        Q_d = x(base + 5);
        
        % Check left boundary: no valid neighbor to the left
        if j == 1 || grid_to_compressed(i, j-1) == 0
            Q_boundary = Q_boundary + Q_l;
        end
        
        % Check right boundary: no valid neighbor to the right
        if j == nx || grid_to_compressed(i, j+1) == 0
            Q_boundary = Q_boundary + Q_r;
        end
        
        % Check top boundary: no valid neighbor above
        if i == 1 || grid_to_compressed(i-1, j) == 0
            Q_boundary = Q_boundary + Q_u;
        end
        
        % Check bottom boundary: no valid neighbor below
        if i == ny || grid_to_compressed(i+1, j) == 0
            Q_boundary = Q_boundary + Q_d;
        end
    end
else
    for c = 1:n_cells
        i = floor((c-1) / nx) + 1;
        j = mod(c-1, nx) + 1;
        
        base = (c-1)*6;
        Q_l = x(base + 2);
        Q_r = x(base + 3);
        Q_u = x(base + 4);
        Q_d = x(base + 5);
        
        % Left boundary
        if j == 1
            Q_boundary = Q_boundary + Q_l;
        end
        % Right boundary
        if j == nx
            Q_boundary = Q_boundary + Q_r;
        end
        % Top boundary
        if i == 1
            Q_boundary = Q_boundary + Q_u;
        end
        % Bottom boundary
        if i == ny
            Q_boundary = Q_boundary + Q_d;
        end
    end
end

end

%% ========================================================================
%  COMPUTE BOUNDARY OUTFLOW FROM CONVERGED STATE
%  ========================================================================

function Q_boundary = compute_boundary_outflow_from_state(x, use_compressed, n_cells, ...
    compressed_to_grid, grid_to_compressed, nx, ny)
% Compute total boundary outflow from the converged state vector x
% This is called AFTER the final update_algebraic_states to ensure
% Q values are consistent with the converged h values

Q_boundary = 0;

if use_compressed
    % For compressed indexing, use compressed_to_grid for fast lookup
    [ny_grid, nx_grid] = size(grid_to_compressed);
    
    for c = 1:n_cells
        % Get (i,j) directly from compressed_to_grid
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
        
        base = (c-1)*6;
        Q_l = x(base + 2);
        Q_r = x(base + 3);
        Q_u = x(base + 4);
        Q_d = x(base + 5);
        
        % Check each direction for boundary (no valid neighbor)
        if j == 1 || grid_to_compressed(i, j-1) == 0
            Q_boundary = Q_boundary + Q_l;
        end
        if j == nx_grid || grid_to_compressed(i, j+1) == 0
            Q_boundary = Q_boundary + Q_r;
        end
        if i == 1 || grid_to_compressed(i-1, j) == 0
            Q_boundary = Q_boundary + Q_u;
        end
        if i == ny_grid || grid_to_compressed(i+1, j) == 0
            Q_boundary = Q_boundary + Q_d;
        end
    end
else
    % For full grid
    for i = 1:ny
        for j = 1:nx
            c = (i-1)*nx + j;
            base = (c-1)*6;
            Q_l = x(base + 2);
            Q_r = x(base + 3);
            Q_u = x(base + 4);
            Q_d = x(base + 5);
            
            if j == 1, Q_boundary = Q_boundary + Q_l; end
            if j == nx, Q_boundary = Q_boundary + Q_r; end
            if i == 1, Q_boundary = Q_boundary + Q_u; end
            if i == ny, Q_boundary = Q_boundary + Q_d; end
        end
    end
end

end

%% Helper function for ternary operator
function result = ternary(condition, true_val, false_val)
if condition
    result = true_val;
else
    result = false_val;
end
end
