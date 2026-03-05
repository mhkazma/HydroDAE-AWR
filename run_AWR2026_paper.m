%% run_AWR2026_paper.m
% Paper reproduction driver for AWR2026:
%   "Evaluating Hydrologic-Hydrodynamic Response to Uncertainty:
%    A Differential-Algebraic Systems Approach"
%
% This script presents the simulation and uncertainty framework.
%   Phase 2: Run DAE simulations + solver comparisons
%   Phase 3: PSE with measurements (+ temporal scheduling)
%   Phase 3b: Monte Carlo validation using independent solvers
%
% USAGE:
%   run_AWR2026_paper            % Run with defaults (vtilted only)
%   WATERSHEDS = {'walnut_gulch'}; run_AWR2026_paper
%   WATERSHEDS = {'vtilted','walnut_gulch'}; run_AWR2026_paper
%
% All outputs go to results/<ws>/AWR2026/.
% Skip logic: existing .mat files are not overwritten unless
% FORCE_RERUN = true.
%
% Author: M. H. Kazma
% Date: February 13, 2026
% Related Paper: AWR2026-HydroStorm

%% ========================================================================
%  CONFIGURATION
% =========================================================================
clear all; close all; clc;

%% --- Global settings ---
% Preserve caller workspace setting for WATERSHEDS if present
if ~exist('WATERSHEDS', 'var')
    WATERSHEDS = {'walnut_gulch'};   % {'vtilted'}, {'walnut_gulch'}, or both
end
if ~exist('FORCE_RERUN', 'var')
    FORCE_RERUN = true;    % If true, all phases will rerun even if output files exist
end

% --- Phase 2: Simulations ---
RUN_DAE_SIMULATION      = true;     % Primary: run_manning_ca_bdf
RUN_SOLVER_COMPARISON   = true;     % CA_Routing + semi-implicit baselines
COMPARISON_SOLVERS      = {'ca_routing', 'semi_implicit'};

% --- Phase 3: PSE Uncertainty ---
RUN_PSE_OPEN_LOOP       = true;    % PSE without measurements (baseline for MC comparison)
RUN_PSE_WITH_MEASUREMENTS = true;    % PSE with measurement matrix (primary)
PSE_METHOD              = 'cell_by_cell';   % 'cell_by_cell' or 'global_sparse'
RUN_UNCERTAINTY_CASES   = false;     % Loop over R, N, KS, PSI, ETA, ALL

% --- Phase 3b: Monte Carlo Validation ---
RUN_MONTE_CARLO         = false;    % Run Monte Carlo validation with independent solvers
MC_N_SAMPLES            = 100;     % Number of MC samples (e.g., 100 for initial validation)
MC_SAMPLING             = 'lhs';    % 'lhs' or 'random'
MC_USE_PARALLEL         = true;
MC_SAVE_TRAJECTORIES    = true;
MC_RANDOM_SEED          = 42;      % Random seed for reproducibility (set to [] for random)
MC_CASES                = {'ALL'};
MC_SOLVERS              = {'ca_routing'};  % or {'ca_routing','semi_implicit'}
MC_LIGHTWEIGHT_SEMI     = true;     

% --- Sensor Configuration ---
SENSOR_STRATEGY         = 'outlet';
SENSOR_UNIFORM_SPACING  = 5;
SENSOR_DEPTH_NOISE_STD  = 0.0005;  % 0.5 mm
SENSOR_FLOW_NOISE_STD   = 1e-4;
SENSOR_INFIL_NOISE_STD  = 0.005;   % 2 mm
SENSOR_OUTLET_FLOW      = true;
SENSOR_OUTLET_NOISE_STD = 5e-3;

% --- Measurement Temporal Scheduling ---
MEASUREMENT_INTERVAL    = 60;       % seconds; [] = continuous

% --- Mesh resolution (vtilted only) ---
VTILTED_DX = 20;  % 5 / 10 / 20 / 40

% --- Comparison solver time steps ---
COMPARISON_CA_DT            = 1.0;
COMPARISON_SEMI_IMPLICIT_DT = 0.5;
COMPARISON_CA_OUTPUT_INTERVAL = 60;

%% ========================================================================
%  SETUP
% =========================================================================
clc; close all;
fprintf('==============================================\n');
fprintf('  AWR2026 Paper\n');
fprintf('==============================================\n');
fprintf('Watersheds: %s\n', strjoin(WATERSHEDS, ', '));
fprintf('Force rerun: %s\n', mat2str(FORCE_RERUN));
fprintf('Date: %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM'));

% Paths
addpath(genpath('src'));

%% ========================================================================
%  MAIN LOOP — iterate over each selected watershed
% =========================================================================

for ws_idx = 1:numel(WATERSHEDS)
    WS = WATERSHEDS{ws_idx};
    fprintf('\n##############################################\n');
    fprintf('  WATERSHED: %s  (%d / %d)\n', upper(WS), ws_idx, numel(WATERSHEDS));
    fprintf('##############################################\n\n');

    % --- Output directory ---
    output_dir = fullfile('results', WS, 'AWR2026');
    if ~isfolder(output_dir), mkdir(output_dir); end
    output_dir_abs = fullfile(pwd, output_dir);

    %% ====================================================================
    %  LOAD PARAMETERS
    % =====================================================================
    fprintf('--- Loading parameters ---\n');
    switch lower(WS)
        case 'vtilted'
            params_vtilted_AWR2026;
            [Z, sim_params] = setup_vtilted_domain(domain.dx, domain.dy, domain.nx, domain.ny);

            sim_params.routing_method = hydraulic.routing_method;
            sim_params.n_base_states = dae.n_base_states;
            sim_params.n_algebraic_states = dae.n_algebraic_states;
            sim_params.n_infiltration_states = dae.n_infiltration_states;
            sim_params.n_states_per_cell = dae.n_states_per_cell;
            sim_params.formulation = dae.formulation;

            if strcmp(infiltration.model, 'green_ampt') && infiltration.enabled
                sim_params = add_infiltration_params(sim_params, ...
                    infiltration.Ks, infiltration.delta_theta, ...
                    infiltration.psi, 0.01);
            end

        case 'walnut_gulch'
            params_walnut_gulch_AWR2026;
            [Z, sim_params] = setup_walnut_gulch_domain(domain.grid_file);

            if strcmp(rainfall.type, 'real_event')
                fl011_row = 70;  fl011_col = 301;
                fl011_ci = sim_params.grid_to_compressed(fl011_row, fl011_col);
                if fl011_ci > 0
                    sim_params.outlet_indices = fl011_ci;
                    sim_params.n_outlet = 1;
                    sim_params.outlet_direction = 'auto';
                    sim_params.slope_outlet = max(0.01, ...
                        sim_params.slope_magnitude(fl011_row, fl011_col));
                    fprintf('  Outlet overridden to FL011 (row %d, col %d, ci=%d)\n', ...
                        fl011_row, fl011_col, fl011_ci);
                else
                    warning('FL011 cell (%d,%d) is not a valid cell — keeping default outlet.', ...
                        fl011_row, fl011_col);
                end
            end

            if isfield(infiltration, 'Ks_factor') && infiltration.Ks_factor ~= 1.0
                sim_params.Ks = sim_params.Ks * infiltration.Ks_factor;
            end
            if isfield(infiltration, 'F_initial')
                sim_params.F_initial = max(infiltration.F_initial, sim_params.F_min);
            end
            if isfield(hydraulic, 'h_smooth')
                sim_params.h_smooth = hydraulic.h_smooth;
            end
            if isfield(hydraulic, 'h_flow_threshold')
                sim_params.h_flow_threshold = hydraulic.h_flow_threshold;
            end

        otherwise
            error('Unknown watershed: %s', WS);
    end

    sim_params.results_dir = output_dir_abs;
    fprintf('  Domain: %d x %d cells, dx=%.0f m\n', ...
        sim_params.grid_size(1), sim_params.grid_size(2), sim_params.dx);

    %% ====================================================================
    %  PHASE 2a: PRIMARY DAE SIMULATION
    % =====================================================================
    result_file = fullfile(output_dir_abs, sprintf('%s_dae_simulation.mat', WS));

    if RUN_DAE_SIMULATION
        fprintf('\n========================================\n');
        fprintf('  PHASE 2a: DAE SIMULATION (%s)\n', upper(dae.solver_method));
        fprintf('========================================\n');

        if ~FORCE_RERUN && isfile(result_file)
            fprintf('  SKIP: Loading existing results from %s\n', result_file);
            loaded = load(result_file, 't_out', 'x_out', 'sim_results', 'sim_params');
            t_out      = loaded.t_out;
            x_out      = loaded.x_out;
            sim_results = loaded.sim_results;
            sim_params  = loaded.sim_params;
            clear loaded;
            % Ensure Q_outlet_ts exists (may be absent in old .mat files)
            sim_results = ensure_Q_outlet_ts(sim_results, x_out, sim_params);
        else
            % --- Prepare rainfall ---
            [R_intensity, duration_s] = prepare_rainfall(rainfall, time, sim_params);

            % --- Dispatch solver ---
            [t_out, x_out, sim_results] = dispatch_dae_solver( ...
                dae, sim_params, R_intensity, duration_s, time, hydraulic);

            fprintf('  Simulation complete: %d time points, %.1f s wallclock\n', ...
                length(t_out), sim_results.computation_time);

            % Water balance verification
            if isfield(sim_results, 'water_balance')
                wb = sim_results.water_balance;
                if isfield(wb, 'rainfall') && wb.rainfall > 0
                    err_pct = 100 * abs(wb.rainfall - wb.infiltration - wb.runoff - wb.ponding) / wb.rainfall;
                    fprintf('  Water balance error: %.3f%%\n', err_pct);
                end
            end

            % Save
            fprintf('  Saving: %s\n', result_file);
            save(result_file, 't_out', 'x_out', 'sim_results', 'sim_params', ...
                'rainfall', 'time', 'hydraulic', 'infiltration', '-v7.3');
        end
    else
        % --- Auto-load previously saved DAE results for downstream phases ---
        if isfile(result_file)
            fprintf('  AUTO-LOAD: DAE results from %s\n', result_file);
            loaded = load(result_file, 't_out', 'x_out', 'sim_results', 'sim_params');
            t_out      = loaded.t_out;
            x_out      = loaded.x_out;
            sim_results = loaded.sim_results;
            sim_params  = loaded.sim_params;
            clear loaded;
            % Ensure Q_outlet_ts exists (may be absent in old .mat files)
            sim_results = ensure_Q_outlet_ts(sim_results, x_out, sim_params);
        end
    end

    %% ====================================================================
    %  PHASE 2b: SOLVER COMPARISON
    % =====================================================================
    if RUN_SOLVER_COMPARISON
        fprintf('\n========================================\n');
        fprintf('  PHASE 2b: SOLVER COMPARISON\n');
        fprintf('========================================\n');

        comp_file = fullfile(output_dir_abs, sprintf('%s_solver_comparison.mat', WS));

        if ~FORCE_RERUN && isfile(comp_file)
            fprintf('  SKIP: %s exists\n', comp_file);
        else
            if ~exist('sim_results', 'var')
                warning('DAE simulation must run first. Skipping comparison.');
            else
                % Prepare rainfall if not already
                if ~exist('R_intensity', 'var')
                    [R_intensity, duration_s] = prepare_rainfall(rainfall, time, sim_params);
                end

                comparison = run_solver_comparison_block( ...
                    sim_params, sim_results, t_out, R_intensity, duration_s, ...
                    dae, COMPARISON_SOLVERS, ...
                    COMPARISON_CA_DT, COMPARISON_SEMI_IMPLICIT_DT, ...
                    COMPARISON_CA_OUTPUT_INTERVAL);

                fprintf('  Saving: %s\n', comp_file);
                save(comp_file, 'comparison', '-v7.3');
            end
        end
    end



    %% ====================================================================
    %  PHASE 3: PSE WITH MEASUREMENTS
    % =====================================================================
    if RUN_PSE_WITH_MEASUREMENTS || RUN_PSE_OPEN_LOOP
        fprintf('\n========================================\n');
        fprintf('  PHASE 3: PSE UNCERTAINTY ANALYSIS\n');
        fprintf('========================================\n');

        if ~exist('sim_results', 'var') || ~exist('x_out', 'var') || isempty(x_out) || ~exist('t_out', 'var')
            warning('DAE simulation with full state trajectory required. Skipping PSE.');
        else
            % --- Configure sensors ---
            sensor_config = configure_sensors(sim_params, ...
                'Strategy', SENSOR_STRATEGY, ...
                'UniformSpacing', SENSOR_UNIFORM_SPACING, ...
                'DepthNoiseStd', SENSOR_DEPTH_NOISE_STD, ...
                'FlowNoiseStd', SENSOR_FLOW_NOISE_STD, ...
                'InfilNoiseStd', SENSOR_INFIL_NOISE_STD, ...
                'OutletFlowSensor', SENSOR_OUTLET_FLOW, ...
                'OutletNoiseStd', SENSOR_OUTLET_NOISE_STD);

            fprintf('  Sensors: %d total (%d depth, %d flow, %d infil)\n', ...
                sensor_config.n_total, sensor_config.n_depth, ...
                sensor_config.n_flow, sensor_config.n_infil);

            % --- Define base PSE config ---
            pse_config = build_pse_config(pse);

            % --- Build temporal measurement schedule ---
            if ~isempty(MEASUREMENT_INTERVAL)
                meas_schedule = build_measurement_schedule(t_out, ...
                    'Interval', MEASUREMENT_INTERVAL);
                pse_config.measurement_schedule = meas_schedule;
                fprintf('  Measurement schedule: every %.0f s (%d/%d timesteps)\n', ...
                    MEASUREMENT_INTERVAL, meas_schedule.n_measurements, meas_schedule.n_total);
            end

            % --- Build uncertainty cases ---
            if RUN_UNCERTAINTY_CASES
                case_list = {'R', 'N', 'KS', 'PSI', 'ETA', 'ALL'};
            else
                case_list = {'ALL'};
            end

            % --- Loop over uncertainty cases ---
            for ci = 1:numel(case_list)
                case_id = case_list{ci};
                pse_case = apply_uncertainty_case(pse_config, case_id);

                % ----- PSE with measurements -----
                if RUN_PSE_WITH_MEASUREMENTS
                    case_file = fullfile(output_dir_abs, ...
                        sprintf('%s_pse_measured_%s.mat', WS, case_id));

                    if ~FORCE_RERUN && isfile(case_file)
                        fprintf('  SKIP: %s (exists)\n', case_file);
                    else
                        fprintf('\n  --- PSE measured, case %s (%s) ---\n', case_id, PSE_METHOD);

                        try
                            switch lower(PSE_METHOD)
                                case 'global_sparse'
                                    pse_meas = run_pse_global_sparse(sim_params, sim_results, ...
                                        t_out, x_out, pse_case, sensor_config);
                                otherwise  % cell_by_cell
                                    pse_meas = run_pse_with_measurements(sim_params, sim_results, ...
                                        t_out, x_out, pse_case, sensor_config);
                            end

                            fprintf('  Done: %.1f s, mean h_std=%.4f mm\n', ...
                                pse_meas.computation_time, ...
                                mean(pse_meas.h_std(:,:,end), 'all') * 1000);

                            save(case_file, 'pse_meas', 'pse_case', 'sensor_config', '-v7.3');
                        catch ME
                            fprintf('  ERROR (PSE %s): %s\n', case_id, ME.message);
                            for kk = 1:min(3, numel(ME.stack))
                                fprintf('    %s L%d\n', ME.stack(kk).name, ME.stack(kk).line);
                            end
                        end
                    end
                end

                % ----- PSE open-loop (no measurements) -----
                if RUN_PSE_OPEN_LOOP
                    ol_file = fullfile(output_dir_abs, ...
                        sprintf('%s_pse_openloop_%s.mat', WS, case_id));

                    if ~FORCE_RERUN && isfile(ol_file)
                        fprintf('  SKIP: %s (exists)\n', ol_file);
                    else
                        fprintf('\n  --- PSE open-loop, case %s ---\n', case_id);
                        try
                            pse_ol = run_pse_dae(sim_params, sim_results, ...
                                t_out, x_out, pse_case);

                            save(ol_file, 'pse_ol', 'pse_case', '-v7.3');
                            fprintf('  Done: %.1f s\n', pse_ol.computation_time);
                        catch ME
                            fprintf('  ERROR (PSE OL %s): %s\n', case_id, ME.message);
                        end
                    end
                end
            end

        end
    end

    %% ====================================================================
    %  PHASE 3b: MONTE CARLO VALIDATION
    % =====================================================================
    if RUN_MONTE_CARLO
        fprintf('\n========================================\n');
        fprintf('  PHASE 3b: MONTE CARLO VALIDATION\n');
        fprintf('========================================\n');

        if ~exist('sim_results', 'var')
            warning('DAE simulation required. Skipping MC.');
        else
            if ~exist('R_intensity', 'var')
                [R_intensity, duration_s] = prepare_rainfall(rainfall, time, sim_params);
            end
            if ~exist('pse_config', 'var')
                pse_config = build_pse_config(pse);
            end

            for mi = 1:numel(MC_CASES)
                mc_case_id = MC_CASES{mi};
                mc_pse_case = apply_uncertainty_case(pse_config, mc_case_id);

                for si = 1:numel(MC_SOLVERS)
                    solver_name = MC_SOLVERS{si};
                    mc_file = fullfile(output_dir_abs, ...
                        sprintf('%s_mc_validation_%s_%s.mat', WS, solver_name, mc_case_id));

                    if ~FORCE_RERUN && isfile(mc_file)
                        fprintf('  SKIP: %s (exists)\n', mc_file);
                    else
                        fprintf('\n  --- MC validation: solver=%s, case=%s, N=%d ---\n', ...
                            solver_name, mc_case_id, MC_N_SAMPLES);
                        try
                            mc_results = run_monte_carlo_validation( ...
                                sim_params, R_intensity, duration_s, ...
                                mc_pse_case, ...
                                'Solver', solver_name, ...
                                'NSamples', MC_N_SAMPLES, ...
                                'Sampling', MC_SAMPLING, ...
                                'UseParallel', MC_USE_PARALLEL, ...
                                'SaveTrajectories', MC_SAVE_TRAJECTORIES, ...
                                'RandomSeed', MC_RANDOM_SEED, ...
                                'LightweightSemi', MC_LIGHTWEIGHT_SEMI, ...
                                'DisplayProgress', true);

                            fprintf('  MC done: %.1f s total, mean σ_h=%.4f mm\n', ...
                                mc_results.computation_time, ...
                                mean(mc_results.h_std_final(:), 'omitnan') * 1000);

                            save(mc_file, 'mc_results', 'mc_pse_case', '-v7.3');
                        catch ME
                            fprintf('  ERROR (MC %s/%s): %s\n', solver_name, mc_case_id, ME.message);
                            for kk = 1:min(3, numel(ME.stack))
                                fprintf('    %s L%d\n', ME.stack(kk).name, ME.stack(kk).line);
                            end
                        end
                    end
                end
            end
        end
    end

    %% Clean up workspace for next watershed iteration
    clear_vars = {'t_out', 'x_out', 'sim_results', 'sim_params', ...
        'R_intensity', 'duration_s', 'comparison', 'sensor_config', ...
        'pse_config', 'pse_meas', 'pse_ol', 'mc_results', ...
        'Z', 'domain', 'hydraulic', 'rainfall', 'infiltration', ...
        'time', 'dae', 'obs', 'estimation', 'pse', 'output', ...
        'viz', 'validation', 'notes', ...
        'initial', 'bc', 'data', 'compute'};
    for cv = 1:numel(clear_vars)
        if exist(clear_vars{cv}, 'var')
            clear(clear_vars{cv});  %#ok<UNRCH>
        end
    end

end % watershed loop

fprintf('\n==============================================\n');
fprintf('  AWR2026 pipeline complete.\n');
fprintf('==============================================\n');


%% ========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function [R_intensity, duration_s] = prepare_rainfall(rainfall, time_cfg, sim_params)
% PREPARE_RAINFALL  Convert rainfall config to solver-ready R_intensity
%
% Replicates the rainfall preparation logic from main_AWR2026.m lines 180–395.

duration_s = time_cfg.end_s;
sim_params.t_end = duration_s;

switch rainfall.type
    case 'uniform_constant'
        R_base = rainfall.intensity / 1000 / 3600;  % mm/hr -> m/s
        rain_dur_s = rainfall.duration * 3600;
        if rain_dur_s < duration_s
            R_intensity = struct('base', R_base, 'duration', rain_dur_s);
        else
            R_intensity = R_base;
        end
        sim_params.rainfall_duration_s = rain_dur_s;

    case 'uniform_schedule'
        schedule_hr = rainfall.schedule;
        schedule_s = schedule_hr;
        schedule_s(:, 1:2) = schedule_hr(:, 1:2) * 3600;
        schedule_s(:, 3) = schedule_hr(:, 3) / 1000 / 3600;

        sim_params.rain_schedule = schedule_s;
        if isfield(rainfall, 'pattern')
            R_intensity = struct('schedule', schedule_s, 'pattern', rainfall.pattern);
        else
            R_intensity = struct('schedule', schedule_s);
        end

        rain_dur_s = 0;
        for rr = 1:size(schedule_s, 1)
            t0 = max(0, schedule_s(rr, 1));
            t1 = min(duration_s, schedule_s(rr, 2));
            if t1 > t0 && schedule_s(rr, 3) > 0
                rain_dur_s = rain_dur_s + (t1 - t0);
            end
        end
        sim_params.rainfall_duration_s = rain_dur_s;

    case 'real_event'
        % Load CSV
        event_table = readtable(rainfall.event_data_file);
        elapsed_min = event_table.ElapsedMin;

        % Sub-window
        event_start = 0;
        event_end = Inf;
        if isfield(rainfall, 'event_start_min'), event_start = rainfall.event_start_min; end
        if isfield(rainfall, 'event_duration_min') && isfinite(rainfall.event_duration_min)
            event_end = event_start + rainfall.event_duration_min;
        end
        keep = (elapsed_min >= event_start) & (elapsed_min < event_end);
        elapsed_min = elapsed_min(keep) - event_start;
        event_table = event_table(keep, :);
        n_intervals = length(elapsed_min);

        use_spatial = ~rainfall.use_watershed_avg;

        if use_spatial
            % IDW spatial rainfall
            gage_nums = rainfall.gage_event_nums;
            n_gages = length(gage_nums);
            gage_data_mmhr = zeros(n_intervals, n_gages);
            for g = 1:n_gages
                col_name = sprintf('Gage_%d_mmhr', gage_nums(g));
                if ismember(col_name, event_table.Properties.VariableNames)
                    gage_data_mmhr(:, g) = event_table.(col_name);
                end
            end

            R_rain_matrix = compute_idw_rainfall( ...
                rainfall.gage_locations, gage_data_mmhr, ...
                sim_params.grid_size, sim_params.dx, sim_params.dy, ...
                sim_params.valid_mask, sim_params.compressed_to_grid, ...
                'Power', rainfall.idw_power, 'DisplayProgress', true);

            R_intensity = struct('rain_matrix', R_rain_matrix, 'dt_rain', 60);
            sim_params.rain_matrix = R_rain_matrix;
            sim_params.dt_rain = 60;
        else
            rain_mmhr = event_table.WatershedAvg_mmhr;
            schedule_s = zeros(n_intervals, 3);
            schedule_s(:, 1) = elapsed_min * 60;
            schedule_s(:, 2) = elapsed_min * 60 + 60;
            schedule_s(:, 3) = rain_mmhr / 1000 / 3600;
            nonzero = schedule_s(:, 3) > 0;
            schedule_s = schedule_s(nonzero, :);
            R_intensity = struct('schedule', schedule_s);
            sim_params.rain_schedule = schedule_s;
        end

        sim_params.rainfall_duration_s = max(elapsed_min) * 60;

    otherwise
        error('Unsupported rainfall type: %s', rainfall.type);
end
end


function [t_out, x_out, sim_results] = dispatch_dae_solver( ...
    dae_cfg, sim_params, R_intensity, duration_s, time_cfg, hydraulic_cfg)
% DISPATCH_DAE_SOLVER  Route to the correct solver based on dae.solver_method
%
% Replicates main_AWR2026.m lines 400–640 (solver dispatch switch).

if strcmp(hydraulic_cfg.routing_method, 'semi_implicit')
    [t_out, sim_results] = run_semi_implicit_simulation(sim_params, R_intensity, duration_s, ...
        'dt', hydraulic_cfg.semi_implicit_dt, 'DisplayProgress', true);
    x_out = [];
    return;
end

% Manning CA with selectable method
solver = dae_cfg.solver_method;
bdf_dt   = safe_field_num(dae_cfg, 'bdf_dt', 0.5);
bdf_ord  = safe_field_num(dae_cfg, 'bdf_order', 1);
bdf_tol  = safe_field_num(dae_cfg, 'bdf_newton_tol', 1e-8);
bdf_iter = safe_field_num(dae_cfg, 'bdf_max_iter', 50);
use_ch   = isfield(dae_cfg, 'use_chordal') && dae_cfg.use_chordal;

switch solver
    case 'bdf'
        fprintf('  Solver: BDF%d fixed-point, dt=%.2f s, tol=%.0e\n', bdf_ord, bdf_dt, bdf_tol);
        [t_out, x_out, sim_results] = run_manning_ca_bdf(sim_params, R_intensity, duration_s, ...
            'dt', bdf_dt, 'BDFOrder', bdf_ord, ...
            'NewtonTol', bdf_tol, 'MaxNewtonIter', bdf_iter, ...
            'DisplayProgress', true);

    case 'bdf_newton'
        fprintf('  Solver: BDF%d Newton [h,F], dt=%.2f s\n', bdf_ord, bdf_dt);
        [t_out, x_out, sim_results] = run_manning_ca_bdf_newton(sim_params, R_intensity, duration_s, ...
            'dt', bdf_dt, 'BDFOrder', bdf_ord, ...
            'NewtonTol', bdf_tol, 'MaxNewtonIter', bdf_iter, ...
            'DisplayProgress', true);

    case 'ode15s'
        ode_args = {'RelTol', time_cfg.RelTol, 'AbsTol', time_cfg.AbsTol, 'DisplayProgress', true};
        if isfield(time_cfg, 'MaxStep'), ode_args = [ode_args, {'MaxStep', time_cfg.MaxStep}]; end
        if isfield(dae_cfg, 'epsilon'), ode_args = [ode_args, {'Epsilon', dae_cfg.epsilon}]; end
        fprintf('  Solver: ode15s semi-explicit DAE\n');
        [t_out, x_out, sim_results] = run_manning_ca_ode15s(sim_params, R_intensity, duration_s, ...
            ode_args{:});

    otherwise
        error('Unknown solver_method: %s', solver);
end
end


function comparison = run_solver_comparison_block( ...
    sim_params, sim_results, t_out, R_intensity, duration_s, ...
    dae_cfg, solvers, ca_dt, semi_dt, ca_output_int)
% RUN_SOLVER_COMPARISON_BLOCK  Run CA_Routing and/or semi-implicit for comparison
%
% Replicates main_AWR2026.m lines 700–1000 (solver comparison block).

ref_label = sprintf('DAE (%s)', dae_cfg.solver_method);
comparison = struct();
comparison.n_solvers = 0;
comparison.ref_label = ref_label;
comparison.ref_t = t_out;
% Use outlet-only (flume) hydrograph for comparison if available
if isfield(sim_results, 'Q_outlet_ts')
    comparison.ref_Q = sim_results.Q_outlet_ts(:);
else
    comparison.ref_Q = sim_results.Q_boundary_ts(:);
end
comparison.ref_h_final = sim_results.h_final;
comparison.ref_wb = sim_results.water_balance;
if isfield(sim_results, 'F_final')
    comparison.ref_F_final = sim_results.F_final;
end
comparison.ref_comp_time = sim_results.computation_time;

% --- CA_Routing ---
if ismember('ca_routing', solvers)
    fprintf('\n  --- CA_Routing (explicit) ---\n');
    try
        ca_results = run_ca_routing_comparison(sim_params, R_intensity, duration_s, ...
            'dt', ca_dt, ...
            'UseInfiltration', true, ...
            'OutputInterval', ca_output_int, ...
            'DisplayProgress', true);

        si = comparison.n_solvers + 1;
        comparison.n_solvers = si;
        comparison.solvers(si).label = 'CA_Routing (Explicit)';
        comparison.solvers(si).t = ca_results.t;
        comparison.solvers(si).Q_outlet = ca_results.Q_outlet;
        comparison.solvers(si).h_final = ca_results.h_final;
        comparison.solvers(si).F_final = ca_results.F_final;
        comparison.solvers(si).water_balance = ca_results.water_balance;
        comparison.solvers(si).computation_time = ca_results.computation_time;
        fprintf('  CA_Routing done: %.2f s\n', ca_results.computation_time);
    catch ME
        fprintf('  ERROR (CA_Routing): %s\n', ME.message);
    end
end

% --- Semi-Implicit ---
if ismember('semi_implicit', solvers)
    fprintf('\n  --- Semi-Implicit (local inertial) ---\n');
    try
        tic;
        [t_semi, semi_res] = run_semi_implicit_simulation(sim_params, R_intensity, duration_s, ...
            'dt', semi_dt, 'DisplayProgress', true);
        semi_comp_time = toc;

        si = comparison.n_solvers + 1;
        comparison.n_solvers = si;
        comparison.solvers(si).label = 'Semi-Implicit (Local Inertial)';
        comparison.solvers(si).t = t_semi;
        if isfield(semi_res, 'Q_outlet_ts')
            comparison.solvers(si).Q_outlet = semi_res.Q_outlet_ts;
        elseif isfield(semi_res, 'Q_boundary_ts')
            comparison.solvers(si).Q_outlet = semi_res.Q_boundary_ts;
        else
            comparison.solvers(si).Q_outlet = [];
        end
        comparison.solvers(si).h_final = semi_res.h(:,:,end);
        if isfield(semi_res, 'F')
            comparison.solvers(si).F_final = semi_res.F(:,:,end);
        else
            comparison.solvers(si).F_final = zeros(size(sim_results.h_final));
        end
        if isfield(semi_res, 'water_balance')
            comparison.solvers(si).water_balance = semi_res.water_balance;
        else
            comparison.solvers(si).water_balance = struct();
        end
        comparison.solvers(si).computation_time = semi_comp_time;
        fprintf('  Semi-Implicit done: %.2f s\n', semi_comp_time);
    catch ME
        fprintf('  ERROR (Semi-Implicit): %s\n', ME.message);
    end
end

% --- Metrics ---
if comparison.n_solvers > 0
    comparison = compute_comparison_metrics(comparison);
end
end


function comparison = compute_comparison_metrics(comparison)
% COMPUTE_COMPARISON_METRICS  NSE, RMSE, R² for each comparison solver

for si = 1:comparison.n_solvers
    if isempty(comparison.solvers(si).Q_outlet), continue; end

    t_ref = comparison.ref_t(:);
    Q_ref = comparison.ref_Q(:);
    t_cmp = comparison.solvers(si).t(:);
    Q_cmp = comparison.solvers(si).Q_outlet(:);

    Q_interp = interp1(t_cmp, Q_cmp, t_ref, 'linear', 'extrap');

    ss_res = sum((Q_ref - Q_interp).^2);
    ss_tot = sum((Q_ref - mean(Q_ref)).^2);
    nse = 1 - ss_res / max(ss_tot, eps);

    rmse_q = sqrt(mean((Q_ref - Q_interp).^2));

    r_val = corrcoef(Q_ref, Q_interp);
    r2 = r_val(1,2)^2;

    comparison.solvers(si).nse = nse;
    comparison.solvers(si).rmse_Q = rmse_q;
    comparison.solvers(si).r2 = r2;

    % Depth RMSE
    h_ref = comparison.ref_h_final(:);
    h_cmp = comparison.solvers(si).h_final(:);
    if numel(h_ref) == numel(h_cmp)
        comparison.solvers(si).rmse_h = sqrt(mean((h_ref - h_cmp).^2));
    end

    fprintf('  %s: NSE=%.4f, RMSE_Q=%.6f m³/s, R²=%.4f\n', ...
        comparison.solvers(si).label, nse, rmse_q, r2);
end

%% AWR2026: Side-by-side water balance comparison
fprintf('\n  ========================================\n');
fprintf('  WATER BALANCE COMPARISON (all solvers)\n');
fprintf('  ========================================\n');
fprintf('  %-28s %10s %10s %10s %10s %8s\n', ...
    'Solver', 'Rain[m³]', 'Infil[m³]', 'Runoff[m³]', 'Pond[m³]', 'Err[%%]');
fprintf('  %-28s %10s %10s %10s %10s %8s\n', ...
    '----------------------------', '----------', '----------', '----------', '----------', '--------');

% Reference (DAE)
wb = comparison.ref_wb;
if isstruct(wb) && isfield(wb, 'rainfall')
    % Compute outflow volume from Q timeseries via trapz
    t_ref = comparison.ref_t(:);
    Q_ref = comparison.ref_Q(:);
    if numel(t_ref) == numel(Q_ref) && numel(t_ref) > 1
        V_out_ref = trapz(t_ref, Q_ref);
    else
        V_out_ref = wb.runoff;
    end
    fprintf('  %-28s %10.1f %10.1f %10.1f %10.1f %7.3f%%\n', ...
        comparison.ref_label, wb.rainfall, wb.infiltration, V_out_ref, wb.ponding, ...
        100*abs(wb.rainfall - wb.infiltration - V_out_ref - wb.ponding)/max(wb.rainfall, 1));
end

for si = 1:comparison.n_solvers
    wb_s = comparison.solvers(si).water_balance;
    if ~isstruct(wb_s) || ~isfield(wb_s, 'rainfall'), continue; end
    
    % Compute outflow volume from Q timeseries
    t_s = comparison.solvers(si).t(:);
    Q_s = comparison.solvers(si).Q_outlet(:);
    if numel(t_s) == numel(Q_s) && numel(t_s) > 1
        V_out_s = trapz(t_s, Q_s);
    elseif isfield(wb_s, 'runoff')
        V_out_s = wb_s.runoff;
    else
        V_out_s = 0;
    end
    
    rain_s = wb_s.rainfall;
    infil_s = wb_s.infiltration;
    if isfield(wb_s, 'ponding')
        pond_s = wb_s.ponding;
    else
        pond_s = rain_s - infil_s - V_out_s;
    end
    
    fprintf('  %-28s %10.1f %10.1f %10.1f %10.1f %7.3f%%\n', ...
        comparison.solvers(si).label, rain_s, infil_s, V_out_s, pond_s, ...
        100*abs(rain_s - infil_s - V_out_s - pond_s)/max(rain_s, 1));
end
fprintf('  ========================================\n');

% Peak Q comparison
fprintf('\n  PEAK DISCHARGE COMPARISON\n');
fprintf('  %-28s %12s %12s\n', 'Solver', 'Q_peak[m³/s]', 't_peak[min]');
fprintf('  %-28s %12s %12s\n', '----------------------------', '------------', '------------');
t_ref = comparison.ref_t(:);
Q_ref = comparison.ref_Q(:);
[Qp, ip] = max(Q_ref);
fprintf('  %-28s %12.4f %12.1f\n', comparison.ref_label, Qp, t_ref(ip)/60);
for si = 1:comparison.n_solvers
    Q_s = comparison.solvers(si).Q_outlet(:);
    t_s = comparison.solvers(si).t(:);
    if isempty(Q_s), continue; end
    [Qp_s, ip_s] = max(Q_s);
    fprintf('  %-28s %12.4f %12.1f\n', comparison.solvers(si).label, Qp_s, t_s(ip_s)/60);
end
fprintf('\n');
end


function pse_config = build_pse_config(pse_struct)
% BUILD_PSE_CONFIG  Create pse_config struct from params file pse struct
pse_config = struct();
pse_config.rainfall_cv      = 0.25;
pse_config.rainfall_sigma_abs = 2.78e-7;
pse_config.manning_n_cv     = 0.20;
pse_config.Ks_cv            = 0.25;
pse_config.psi_cv           = 0.10;
pse_config.eta_cv           = 0.12;
pse_config.confidence_level = 0.95;
pse_config.debug            = true;
pse_config.stride           = 1;

% Dry cell fast-path configuration
pse_config.dry_cell_threshold   = 1e-8;   % h below this triggers fast path [m]
pse_config.dry_cell_decay_rain  = 0.5;    % variance decay factor when raining
pse_config.dry_cell_decay_norain = 0.1;   % variance decay factor when dry

% Stride interpolation method: 'linear', 'pchip', or 'hold'
pse_config.stride_interp = 'linear';

if exist('pse_struct', 'var') && ~isempty(pse_struct) && isstruct(pse_struct)
    if isfield(pse_struct, 'distribution_type')
        pse_config.distribution_type = pse_struct.distribution_type;
    else
        pse_config.distribution_type = 'gaussian';
    end
    if isfield(pse_struct, 'param_distributions')
        pse_config.param_distributions = pse_struct.param_distributions;
    end
    if isfield(pse_struct, 'param_bounds')
        pse_config.param_bounds = pse_struct.param_bounds;
    end
    if isfield(pse_struct, 'min_depth_for_gaussian')
        pse_config.min_depth_for_gaussian = pse_struct.min_depth_for_gaussian;
    end
    if isfield(pse_struct, 'truncation_lower')
        pse_config.truncation_lower = pse_struct.truncation_lower;
    end
else
    pse_config.distribution_type = 'gaussian';
end
end


function pse_case = apply_uncertainty_case(pse_config, case_id)
% APPLY_UNCERTAINTY_CASE  Zero out CVs except the target source
%
% Cases: R, N, KS, PSI, ETA, ALL

pse_case = pse_config;  % copy

switch upper(case_id)
    case 'R'
        pse_case.manning_n_cv = 0;
        pse_case.Ks_cv = 0;
        pse_case.psi_cv = 0;
        pse_case.eta_cv = 0;
    case 'N'
        pse_case.rainfall_cv = 0;
        pse_case.rainfall_sigma_abs = 0;
        pse_case.Ks_cv = 0;
        pse_case.psi_cv = 0;
        pse_case.eta_cv = 0;
    case 'KS'
        pse_case.rainfall_cv = 0;
        pse_case.rainfall_sigma_abs = 0;
        pse_case.manning_n_cv = 0;
        pse_case.psi_cv = 0;
        pse_case.eta_cv = 0;
    case 'PSI'
        pse_case.rainfall_cv = 0;
        pse_case.rainfall_sigma_abs = 0;
        pse_case.manning_n_cv = 0;
        pse_case.Ks_cv = 0;
        pse_case.eta_cv = 0;
    case 'ETA'
        pse_case.rainfall_cv = 0;
        pse_case.rainfall_sigma_abs = 0;
        pse_case.manning_n_cv = 0;
        pse_case.Ks_cv = 0;
        pse_case.psi_cv = 0;
    case 'ALL'
        % Keep all CVs as-is
    otherwise
        warning('Unknown uncertainty case: %s. Using ALL.', case_id);
end
end


function val = safe_field_num(s, field, default)
% SAFE_FIELD_NUM  Get numeric field with default
if isfield(s, field)
    val = s.(field);
else
    val = default;
end
end


function sim_results = ensure_Q_outlet_ts(sim_results, x_out, sim_params)
% ENSURE_Q_OUTLET_TS  Compute Q_outlet_ts from x_out if not already present.
% For old .mat files that only have Q_boundary_ts, extract outlet-only Q
% from the DAE state trajectory (6 states per cell: h, Q_l, Q_r, Q_u, Q_d, F).

if isfield(sim_results, 'Q_outlet_ts')
    return;  % Already computed
end

if ~isfield(sim_params, 'outlet_indices') || isempty(sim_params.outlet_indices)
    % No specific outlets — fall back to Q_boundary_ts
    if isfield(sim_results, 'Q_boundary_ts')
        sim_results.Q_outlet_ts = sim_results.Q_boundary_ts;
    end
    return;
end

outlet_indices = sim_params.outlet_indices;
use_compressed = isfield(sim_params, 'use_compressed_indexing') && sim_params.use_compressed_indexing;
n_times = size(x_out, 1);

if use_compressed
    ctg = sim_params.compressed_to_grid;
    gtc = sim_params.grid_to_compressed;
    nx = sim_params.grid_size(2);
    ny = sim_params.grid_size(1);
else
    nx = sim_params.grid_size(2);
    ny = sim_params.grid_size(1);
end

Q_outlet_ts = zeros(n_times, 1);
for ti = 1:n_times
    x = x_out(ti, :)';
    Q_outlet = 0;
    for k = 1:length(outlet_indices)
        c = outlet_indices(k);
        base = (c-1)*6;
        if use_compressed
            i = ctg(c, 1);
            j = ctg(c, 2);
            if j == 1 || gtc(i, j-1) == 0
                Q_outlet = Q_outlet + x(base + 2);
            end
            if j == nx || gtc(i, j+1) == 0
                Q_outlet = Q_outlet + x(base + 3);
            end
            if i == 1 || gtc(i-1, j) == 0
                Q_outlet = Q_outlet + x(base + 4);
            end
            if i == ny || gtc(i+1, j) == 0
                Q_outlet = Q_outlet + x(base + 5);
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
    Q_outlet_ts(ti) = Q_outlet;
end

sim_results.Q_outlet_ts = Q_outlet_ts;
fprintf('  Computed Q_outlet_ts from x_out (%d outlet cells, %d timesteps)\n', ...
    length(outlet_indices), n_times);
end
