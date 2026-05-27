% AWR2026_run_mc_celltrack.m
% =========================================================================
% Author: M. H. Kazma
% Date: May 2026
% =========================================================================

clear; clc;
fprintf('\n======================================================\n');
fprintf('  AWR2026 R1 C1-R2: Tracked-cell MC validation\n');
fprintf('======================================================\n');

% ---- Path setup (mirrors run_AWR2026_paper.m) ----
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_dir);
cd(repo_root);
addpath(genpath(fullfile(repo_root, 'src')));
addpath(repo_root);
% HydroPol2D CA_Routing dependency (external; mirrors interactive MATLAB userpath)
hp2d_root = fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'HydroPol2D');
if isfolder(hp2d_root)
    addpath(genpath(hp2d_root));
end

% ---- Shared MC settings (must match run_AWR2026_paper.m) ----
MC_N_SAMPLES        = 100;
MC_SAMPLING         = 'lhs';
MC_USE_PARALLEL     = true;
MC_SAVE_TRAJECTORIES = true;
MC_RANDOM_SEED      = 42;
MC_LIGHTWEIGHT_SEMI = false;

WATERSHEDS = {'vtilted', 'walnut_gulch'};

for ws_idx = 1:numel(WATERSHEDS)
    WS = WATERSHEDS{ws_idx};
    fprintf('\n##############################################\n');
    fprintf('  WATERSHED: %s  (%d / %d)\n', upper(WS), ws_idx, numel(WATERSHEDS));
    fprintf('##############################################\n\n');

    output_dir = fullfile('results', WS, 'AWR2026');
    if ~isfolder(output_dir), mkdir(output_dir); end
    output_dir_abs = fullfile(pwd, output_dir);

    mc_file = fullfile(output_dir_abs, ...
        sprintf('%s_mc_celltrack_ca_routing_ALL.mat', WS));
    if isfile(mc_file)
        fprintf('  EXISTS: %s\n  Skipping (delete file to redo).\n', mc_file);
        continue;
    end

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

            % V-Tilted: obs.locations stored as [col, row]; swap.
            track_cells_local = [obs.locations(:, 2), obs.locations(:, 1)];
            track_cells_local = max(track_cells_local, 1);
            track_cells_local(:, 1) = min(track_cells_local(:, 1), domain.ny);
            track_cells_local(:, 2) = min(track_cells_local(:, 2), domain.nx);
            track_labels_local = arrayfun(@(k) sprintf('obs_%d', k), ...
                1:size(track_cells_local, 1), 'UniformOutput', false);

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

            % WGEW: 5 gauged flumes + 2 ungauged interior cells.
            % Hardcode flume (row,col) on the 242x521 grid (obs network is
            % only auto-configured when rainfall.type == 'real_event').
            gauged = [100,   1;   % FL001
                       87, 104;   % FL002
                      108, 185;   % FL006
                      121, 243;   % FL009
                       70, 301];  % FL011 (outlet)
            gauge_names_local = {'FL001','FL002','FL006','FL009','FL011'};
            ny_wg = sim_params.grid_size(1);
            nx_wg = sim_params.grid_size(2);

            ungauged = [round(0.55 * ny_wg), round(0.40 * nx_wg);
                        round(0.65 * ny_wg), round(0.75 * nx_wg)];
            ungauged = max(ungauged, 1);
            ungauged(:, 1) = min(ungauged(:, 1), ny_wg);
            ungauged(:, 2) = min(ungauged(:, 2), nx_wg);

            if isfield(sim_params, 'valid_mask')
                vm = sim_params.valid_mask;
                for k = 1:size(ungauged, 1)
                    r = ungauged(k, 1);  c = ungauged(k, 2);
                    if ~vm(r, c)
                        found = false;
                        for radius = 1:max(ny_wg, nx_wg)
                            for dr = -radius:radius
                                for dc = -radius:radius
                                    rr = min(max(r + dr, 1), ny_wg);
                                    cc = min(max(c + dc, 1), nx_wg);
                                    if vm(rr, cc)
                                        ungauged(k, :) = [rr, cc];
                                        found = true;  break;
                                    end
                                end
                                if found, break; end
                            end
                            if found, break; end
                        end
                    end
                end
            end

            track_cells_local = [gauged; ungauged];
            track_labels_local = [gauge_names_local, {'ungauged_A', 'ungauged_B'}];

        otherwise
            error('Unknown watershed: %s', WS);
    end

    sim_params.results_dir = output_dir_abs;
    fprintf('  Domain: %d x %d cells, dx=%.1f m\n', ...
        sim_params.grid_size(1), sim_params.grid_size(2), sim_params.dx);

    fprintf('  Track cells (row, col):\n');
    for k = 1:size(track_cells_local, 1)
        fprintf('    %-12s : (%d, %d)\n', track_labels_local{k}, ...
            track_cells_local(k, 1), track_cells_local(k, 2));
    end

    % ---- Rainfall + PSE config ----
    [R_intensity, duration_s] = prepare_rainfall(rainfall, time, sim_params);
    pse_config = build_pse_config(pse);

    % Match manuscript MC: case_id = 'ALL' (all uncertainty sources active)
    mc_pse_case = apply_uncertainty_case(pse_config, 'ALL');

    if strcmpi(WS, 'walnut_gulch')
        mc_pse_case.rainfall_cv     = 0.20;
        mc_pse_case.rainfall_sigma_abs = 0;
        mc_pse_case.manning_n_cv    = 0.20;
        mc_pse_case.Ks_cv           = 0.20;
        mc_pse_case.psi_cv          = 0.20;
        mc_pse_case.eta_cv          = 0.20;
        fprintf('  [WGEW override] All 5 CVs forced to 20%% to match manuscript MC.\n');
    end

    fprintf('\n  Running MC validation: N=%d, seed=%d, n_track=%d\n', ...
        MC_N_SAMPLES, MC_RANDOM_SEED, size(track_cells_local, 1));
    t0 = tic;
    mc_results = run_monte_carlo_validation( ...
        sim_params, R_intensity, duration_s, mc_pse_case, ...
        'Solver', 'ca_routing', ...
        'NSamples', MC_N_SAMPLES, ...
        'Sampling', MC_SAMPLING, ...
        'UseParallel', MC_USE_PARALLEL, ...
        'SaveTrajectories', MC_SAVE_TRAJECTORIES, ...
        'RandomSeed', MC_RANDOM_SEED, ...
        'LightweightSemi', MC_LIGHTWEIGHT_SEMI, ...
        'TrackCells', track_cells_local, ...
        'DisplayProgress', true);
    fprintf('  Done: %.1f s, %d/%d successful\n', ...
        toc(t0), mc_results.n_successful, MC_N_SAMPLES);

    track_cells  = track_cells_local;
    track_labels = track_labels_local;
    save(mc_file, 'mc_results', 'track_cells', 'track_labels', ...
        'mc_pse_case', '-v7.3');
    fprintf('  Saved: %s\n', mc_file);

    clear_vars = {'sim_params', 'R_intensity', 'duration_s', 'pse_config', ...
        'mc_pse_case', 'mc_results', 'Z', 'domain', 'hydraulic', ...
        'rainfall', 'infiltration', 'time', 'dae', 'obs', 'estimation', ...
        'pse', 'output', 'viz', 'validation', 'notes', 'initial', 'bc', ...
        'track_cells', 'track_labels', 'track_cells_local', 'track_labels_local'};
    for cv = 1:numel(clear_vars)
        if exist(clear_vars{cv}, 'var')
            clear(clear_vars{cv});
        end
    end
end

fprintf('\n======================================================\n');
fprintf('  Tracked-cell MC complete.\n');
fprintf('======================================================\n');


% =========================================================================
% LOCAL FUNCTIONS (copied verbatim from run_AWR2026_paper.m to keep the
% wrapper self-contained without re-running the full driver).
% =========================================================================

function [R_intensity, duration_s] = prepare_rainfall(rainfall, time_cfg, sim_params)
duration_s = time_cfg.end_s;
sim_params.t_end = duration_s; %#ok<NASGU>

switch rainfall.type
    case 'uniform_constant'
        R_base = rainfall.intensity / 1000 / 3600;
        rain_dur_s = rainfall.duration * 3600;
        if rain_dur_s < duration_s
            R_intensity = struct('base', R_base, 'duration', rain_dur_s);
        else
            R_intensity = R_base;
        end

    case 'uniform_schedule'
        schedule_hr = rainfall.schedule;
        schedule_s = schedule_hr;
        schedule_s(:, 1:2) = schedule_hr(:, 1:2) * 3600;
        schedule_s(:, 3) = schedule_hr(:, 3) / 1000 / 3600;
        if isfield(rainfall, 'pattern')
            R_intensity = struct('schedule', schedule_s, 'pattern', rainfall.pattern);
        else
            R_intensity = struct('schedule', schedule_s);
        end

    case 'real_event'
        event_table = readtable(rainfall.event_data_file);
        elapsed_min = event_table.ElapsedMin;

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
        else
            rain_mmhr = event_table.WatershedAvg_mmhr;
            schedule_s = zeros(n_intervals, 3);
            schedule_s(:, 1) = elapsed_min * 60;
            schedule_s(:, 2) = elapsed_min * 60 + 60;
            schedule_s(:, 3) = rain_mmhr / 1000 / 3600;
            nonzero = schedule_s(:, 3) > 0;
            schedule_s = schedule_s(nonzero, :);
            R_intensity = struct('schedule', schedule_s);
        end

    otherwise
        error('Unsupported rainfall type: %s', rainfall.type);
end
end


function pse_config = build_pse_config(pse_struct)
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
pse_config.dry_cell_threshold    = 1e-8;
pse_config.dry_cell_decay_rain   = 0.5;
pse_config.dry_cell_decay_norain = 0.1;
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
pse_case = pse_config;
switch upper(case_id)
    case 'R'
        pse_case.manning_n_cv = 0; pse_case.Ks_cv = 0;
        pse_case.psi_cv = 0; pse_case.eta_cv = 0;
    case 'N'
        pse_case.rainfall_cv = 0; pse_case.rainfall_sigma_abs = 0;
        pse_case.Ks_cv = 0; pse_case.psi_cv = 0; pse_case.eta_cv = 0;
    case 'KS'
        pse_case.rainfall_cv = 0; pse_case.rainfall_sigma_abs = 0;
        pse_case.manning_n_cv = 0; pse_case.psi_cv = 0; pse_case.eta_cv = 0;
    case 'PSI'
        pse_case.rainfall_cv = 0; pse_case.rainfall_sigma_abs = 0;
        pse_case.manning_n_cv = 0; pse_case.Ks_cv = 0; pse_case.eta_cv = 0;
    case 'ETA'
        pse_case.rainfall_cv = 0; pse_case.rainfall_sigma_abs = 0;
        pse_case.manning_n_cv = 0; pse_case.Ks_cv = 0; pse_case.psi_cv = 0;
    case 'ALL'
        % keep all
    otherwise
        warning('Unknown uncertainty case: %s. Using ALL.', case_id);
end
end

