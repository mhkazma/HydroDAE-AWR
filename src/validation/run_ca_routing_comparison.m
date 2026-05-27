function results = run_ca_routing_comparison(sim_params, R_input, duration, varargin)
% RUN_CA_ROUTING_COMPARISON  Run original HydroPol2D CA_Routing for comparison
%
% Runs the explicit CA_Routing solver on the same
% domain and forcing as the DAE solver for verification/comparison.
% Includes optional Green-Ampt infiltration implemented inline (same GA
% formula as the BDF solver for fair comparison).
%
% Inputs:
%   sim_params - Domain parameters struct (from setup_vtilted_domain, etc.)
%   R_input    - Rainfall forcing:
%                  scalar [m/s]: constant rainfall
%                  struct with .base [m/s] and .duration [s]: time-limited
%                  struct with .schedule [Nx3: t_start_s, t_end_s, R_m/s]
%                  struct with .rain_matrix and .dt_rain: spatiotemporal
%   duration   - Simulation duration [s]
%
% Optional Name-Value:
%   'dt'                - Time step [s] (default: 1.0)
%   'UseInfiltration'   - Enable Green-Ampt (default: true if params.Ks exists)
%   'DisplayProgress'   - Print progress (default: true)
%   'OutputInterval'    - Save output every N steps (default: 60)
%
% Outputs:
%   results.t                  - Time vector [s]
%   results.Q_outlet           - Outlet discharge [m^3/s]
%   results.h_final            - Final depth map [m]
%   results.F_final            - Final cumulative infiltration [m]
%   results.water_balance      - Struct: .rainfall, .infiltration, .ponding, .runoff [m^3]
%   results.computation_time   - Wall clock [s]
%   results.solver_name        - 'CA_Routing'
%
% Notes:
%   - CA_Routing works internally in mm (depths) and minutes (timestep).
%   - outlet_flow from CA_Routing is in mm/hr; converted to m^3/s here.
%   - Entire bottom row is marked as outlet (matching BDF solver boundary).
%   - Green-Ampt: f = Ks * (1 + eta * (h + psi_f) / max(F, F_min))
%
% Author: M. H. Kazma
% Date: February 9, 2026

%% Parse inputs
p = inputParser;
addRequired(p, 'sim_params');
addRequired(p, 'R_input');
addRequired(p, 'duration', @isnumeric);
addParameter(p, 'dt', 1.0, @isnumeric);
addParameter(p, 'UseInfiltration', [], @islogical);
addParameter(p, 'DisplayProgress', true, @islogical);
addParameter(p, 'OutputInterval', 60, @isnumeric);
addParameter(p, 'TrackCells', [], @isnumeric);  % [n_track x 2] (row, col) pairs to record per-step depth
parse(p, sim_params, R_input, duration, varargin{:});

dt_s = p.Results.dt;
display_progress = p.Results.DisplayProgress;
output_interval = p.Results.OutputInterval;
track_cells = p.Results.TrackCells;
if ~isempty(track_cells) && size(track_cells, 2) ~= 2
    error('run_ca_routing_comparison: TrackCells must be [n_track x 2] (row, col).');
end
n_track = size(track_cells, 1);

% Infiltration: default true if Ks is available
use_infiltration = p.Results.UseInfiltration;
if isempty(use_infiltration)
    use_infiltration = isfield(sim_params, 'Ks');
end

%% Extract domain parameters
ny = sim_params.grid_size(1);
nx = sim_params.grid_size(2);
dx = sim_params.dx;
Z = sim_params.Z;
n_manning = sim_params.n_manning;
cell_area = dx * dx;     % CA_Routing assumes square cells
n_cells = ny * nx;

% Initial abstraction (convert m -> mm for CA_Routing)
if isfield(sim_params, 'h_0')
    h_0_mm = sim_params.h_0 * 1000;
else
    h_0_mm = 0;
end
h_0_cell = h_0_mm * ones(ny, nx);

% Outlet configuration
slope_outlet = sim_params.slope_outlet;
outlet_type = 1;  % Normal depth (default)
if isfield(sim_params, 'outlet_type')
    outlet_type = sim_params.outlet_type;
end

% No-data mask (needed before outlet_check expansion in Fix 5)
idx_nan = false(ny, nx);
if isfield(sim_params, 'valid_mask')
    idx_nan = ~sim_params.valid_mask;
end

% Build outlet_check matrix (mark outlet cells)
outlet_check = zeros(ny, nx);
use_compressed = isfield(sim_params, 'use_compressed_indexing') && sim_params.use_compressed_indexing;
if isfield(sim_params, 'outlet_direction') && strcmp(sim_params.outlet_direction, 'down')
    % vtilted-style: entire bottom row as outlet
    outlet_check(ny, :) = 1;
elseif isfield(sim_params, 'outlet_indices')
    if use_compressed
        % Compressed indexing (e.g., walnut_gulch): convert via lookup table
        ctg = sim_params.compressed_to_grid;
        for k = 1:length(sim_params.outlet_indices)
            ci = sim_params.outlet_indices(k);
            outlet_check(ctg(ci, 1), ctg(ci, 2)) = 1;
        end
    else
        % Full-grid linear indices (e.g., vtilted)
        [oi, oj] = ind2sub([ny, nx], sim_params.outlet_indices);
        for k = 1:length(oi)
            outlet_check(oi(k), oj(k)) = 1;
        end
    end
else
    outlet_check(ny, :) = 1;  % Default: bottom row
end

outlet_measurement = outlet_check;  % Save original (flume/outlet cells only)
if isfield(sim_params, 'is_boundary')
    is_bnd = sim_params.is_boundary;  % (ny, nx) logical or double
    if size(is_bnd, 1) ~= ny
        is_bnd = is_bnd';  % Handle h5py transpose
    end
    n_outlet_orig = sum(outlet_check(:));
    outlet_check = double(is_bnd > 0);
    % Zero out invalid (outside-watershed) cells from outlet_check
    outlet_check(idx_nan) = 0;
    n_outlet_new = sum(outlet_check(:));
    if display_progress
        fprintf('  FIX 5: outlet_check expanded from %d → %d boundary cells (drainage)\n', ...
            n_outlet_orig, n_outlet_new);
        fprintf('         outlet_measurement kept at %d cells (Q_out tracking)\n', ...
            sum(outlet_measurement(:)));
    end
end
% ---- END FIX 5 ----

% Row/col of outlet cells (needed by CA_Routing for outlet_type == 2)
[row_outlet, col_outlet] = find(outlet_check);

% Velocity constraint flag (0 = no critical velocity limit = HydroPol2D a)
flag_critical = 0;

% Reservoir params (unused - pass empty)
reservoir_x = []; reservoir_y = [];
k1 = []; h1 = []; k2 = []; k3 = []; h2 = []; k4 = [];
yds1 = []; xds1 = []; yds2 = []; xds2 = [];
flag_reservoir = 0;

%% Parse rainfall forcing
if isstruct(R_input)
    if isfield(R_input, 'schedule')
        rain_type = 'schedule';
        schedule = R_input.schedule;  % [t_start_s, t_end_s, R_m_per_s]
    elseif isfield(R_input, 'base') && isfield(R_input, 'duration')
        rain_type = 'time_limited';
        R_base = R_input.base;
        rain_duration = R_input.duration;
    elseif isfield(R_input, 'rain_matrix')
        rain_type = 'matrix';
        rain_matrix = R_input.rain_matrix;
        dt_rain = R_input.dt_rain;
    else
        error('run_ca_routing_comparison: Unrecognized R_input struct format.');
    end
else
    rain_type = 'constant';
    R_const = R_input;  % m/s
end

%% Infiltration parameters
if use_infiltration
    if ~isfield(sim_params, 'Ks')
        error('run_ca_routing_comparison: Infiltration requested but sim_params.Ks not found.');
    end
    % Expand scalars to matrices if needed
    if isscalar(sim_params.Ks)
        Ks = sim_params.Ks * ones(ny, nx);
    else
        Ks = sim_params.Ks;
    end
    if isscalar(sim_params.eta)
        eta = sim_params.eta * ones(ny, nx);
    else
        eta = sim_params.eta;
    end
    if isscalar(sim_params.psi_f)
        psi_f = sim_params.psi_f * ones(ny, nx);
    else
        psi_f = sim_params.psi_f;
    end
    F_min = sim_params.F_min;
    F_initial = F_min;
    if isfield(sim_params, 'F_initial')
        F_initial = sim_params.F_initial;
    end
end

%% Initialize state
d_tot = zeros(ny, nx);              % Depth [mm] (CA_Routing convention)
I_tot_end_cell = zeros(ny, nx);     % Previous outflow volume [m^3]

if use_infiltration
    F_cum = F_initial * ones(ny, nx);  % Cumulative infiltration [m]
end

%% Time stepping
time_step_min = dt_s / 60;  % CA_Routing expects time_step in minutes
n_steps = ceil(duration / dt_s);
n_output = ceil(n_steps / output_interval) + 2;

t_save = zeros(n_output, 1);
Q_outlet_save = zeros(n_output, 1);
if n_track > 0
    h_track_save = zeros(n_track, n_output);  % depth [m] at tracked cells
else
    h_track_save = [];
end
save_idx = 1;
t_save(1) = 0;
Q_outlet_save(1) = 0;

total_rainfall_vol = 0;      % m^3
total_infiltration_vol = 0;  % m^3
total_outflow_vol = 0;       % m^3

if display_progress
    fprintf('  Grid: %d x %d = %d cells\n', ny, nx, n_cells);
    fprintf('  Time step: %.2f s (%.4f min)\n', dt_s, time_step_min);
    fprintf('  Total steps: %d\n', n_steps);
    fprintf('  Duration: %.0f s (%.1f min)\n', duration, duration/60);
    if use_infiltration
        if isfield(sim_params, 'valid_mask')
            fprintf('  Infiltration: Green-Ampt (Ks=%.2f mm/h, mean over valid cells)\n', ...
                nanmean(Ks(sim_params.valid_mask)) * 3600 * 1000);
        else
            fprintf('  Infiltration: Green-Ampt (Ks=%.2f mm/h)\n', Ks(1,1)*3600*1000);
        end
    else
        fprintf('  Infiltration: disabled\n');
    end
    fprintf('\n  Running CA_Routing time loop...\n');
end

%% Main time loop
tic;
for step = 1:n_steps
    t = (step - 1) * dt_s;
    t_new = min(step * dt_s, duration);
    dt_actual = t_new - t;
    
    %% Compute rainfall at current time [m/s]
    switch rain_type
        case 'constant'
            R_now = R_const;
        case 'time_limited'
            if t < rain_duration
                R_now = R_base;
            else
                R_now = 0;
            end
        case 'schedule'
            R_now = 0;
            for rr = 1:size(schedule, 1)
                if t >= schedule(rr, 1) && t < schedule(rr, 2)
                    R_now = schedule(rr, 3);
                    break;
                end
            end
        case 'matrix'
            time_idx = min(floor(t / dt_rain) + 1, size(rain_matrix, 2));
            if use_compressed && isfield(sim_params, 'compressed_to_grid')
                % Expand compressed rain vector to full grid (per-cell spatial rainfall)
                R_now_field = zeros(ny, nx);
                ctg_rain = sim_params.compressed_to_grid;
                for rc = 1:size(ctg_rain, 1)
                    R_now_field(ctg_rain(rc, 1), ctg_rain(rc, 2)) = rain_matrix(rc, time_idx);
                end
            else
                R_now_field = mean(rain_matrix(:, time_idx)) * ones(ny, nx);
            end
            R_now = mean(rain_matrix(:, time_idx));  % Scalar for diagnostics
    end
    
    % Add rainfall to depth (m/s * s = m -> mm)
    if exist('R_now_field', 'var') && strcmp(rain_type, 'matrix')
        rainfall_mm_field = R_now_field * dt_actual * 1000;
        d_tot = d_tot + rainfall_mm_field;
        total_rainfall_vol = total_rainfall_vol + sum(R_now_field(:)) * dt_actual * cell_area;
    else
        rainfall_mm = R_now * dt_actual * 1000;
        if any(idx_nan(:))
            d_tot(~idx_nan) = d_tot(~idx_nan) + rainfall_mm;
        else
            d_tot = d_tot + rainfall_mm;
        end
        
        if isfield(sim_params, 'valid_mask')
            total_rainfall_vol = total_rainfall_vol + R_now * dt_actual * sum(sim_params.valid_mask(:)) * cell_area;
        else
            total_rainfall_vol = total_rainfall_vol + R_now * dt_actual * n_cells * cell_area;
        end
    end
    
    %% Green-Ampt infiltration (before routing)
    if use_infiltration
        h_m = d_tot / 1000;  % Current depth [m]
        
        % Infiltration capacity [m/s]: f = Ks * (1 + eta*(h + psi_f) / max(F, F_min))
        f_cap = Ks .* (1 + eta .* (h_m + psi_f) ./ max(F_cum, F_min));
        
        % Actual infiltration: limited by available water
        f_actual = min(f_cap, h_m / dt_actual);
        
        % Only infiltrate where there is water
        f_actual(h_m < 1e-8) = 0;
        
        % ORIGINAL (HydroPol2D):
        f_actual(idx_nan) = 0;
        
        % Update cumulative infiltration and depth
        dF = f_actual * dt_actual;  % [m]
        F_cum = F_cum + dF;
        d_tot = d_tot - dF * 1000;  % Remove infiltrated depth [mm]
        d_tot = max(d_tot, 0);
        
        % ORIGINAL (HydroPol2D):
        %   total_infiltration_vol = total_infiltration_vol + sum(dF(:)) * cell_area;
        % This summed infiltration over ALL cells including invalid ones.
        %
        % FIX [AWR2026]: Sum infiltration volume only over valid cells for
        % consistency with rainfall volume tracking.
        if any(idx_nan(:))
            total_infiltration_vol = total_infiltration_vol + sum(dF(~idx_nan)) * cell_area;
        else
            total_infiltration_vol = total_infiltration_vol + sum(dF(:)) * cell_area;
        end
    end
    
    %% CA_Routing (overland flow routing)
    [qout_left, qout_right, qout_up, qout_down, outlet_flow, d_t, I_tot_end_cell] = ...
        CA_Routing(reservoir_x, reservoir_y, k1, h1, k2, k3, h2, k4, ...
                   yds1, xds1, yds2, xds2, flag_reservoir, ...
                   Z, d_tot, n_manning, cell_area, time_step_min, ...
                   h_0_cell, dx, I_tot_end_cell, ...
                   outlet_check, outlet_type, slope_outlet, ...
                   row_outlet, col_outlet, idx_nan, flag_critical);
    
    %% Outflows become Inflows (critical redistribution step)
    % Neighbor's outflow in opposite direction becomes this cell's inflow
    % qout is in mm/hr; multiply by (time_step_min / 60) to get mm
    qin_left  = [zeros(ny, 1), qout_right(:, 1:(nx-1))];   % right-flow from left neighbor
    qin_right = [qout_left(:, (2:nx)), zeros(ny, 1)];       % left-flow from right neighbor
    qin_up    = [zeros(1, nx); qout_down(1:(ny-1), :)];     % down-flow from upper neighbor
    qin_down  = [qout_up(2:ny, :); zeros(1, nx)];           % up-flow from lower neighbor
    qin_total = qin_left + qin_right + qin_up + qin_down;   % mm/hr
    
    % Clean NaN/Inf
    qin_total(isnan(qin_total) | isinf(qin_total)) = 0;
    
    % Update depth: d_t already has outflow subtracted; add neighbor inflow
    d_t = d_t + qin_total * (time_step_min / 60);  % mm/hr * hr = mm
    d_t = max(d_t, 0);  % Prevent negative depths
    
    % Update depth
    d_tot = d_t;
    
    % Extract outlet discharge [m^3/s]
    % outlet_flow is in mm/hr per cell; convert: mm/hr * m^2 / (1000 * 3600) = m^3/s
    % FIX [AWR2026]: Use outlet_measurement (flume cells) for Q_out,
    % matching DAE Q_outlet_ts which tracks designated outlet cells only.
    Q_out_m3s = sum(outlet_flow(:) .* outlet_measurement(:)) * cell_area / (1000 * 3600);
    total_outflow_vol = total_outflow_vol + Q_out_m3s * dt_actual;
    
    %% Save output at intervals
    if mod(step, output_interval) == 0 || step == n_steps
        save_idx = save_idx + 1;
        if save_idx > n_output
            t_save = [t_save; zeros(200, 1)];
            Q_outlet_save = [Q_outlet_save; zeros(200, 1)];
            if n_track > 0
                h_track_save = [h_track_save, zeros(n_track, 200)];
            end
            n_output = n_output + 200;
        end
        t_save(save_idx) = t_new;
        Q_outlet_save(save_idx) = Q_out_m3s;
        if n_track > 0
            for kk = 1:n_track
                h_track_save(kk, save_idx) = d_tot(track_cells(kk,1), track_cells(kk,2)) / 1000;
            end
        end
    end
    
    %% Progress updates
    if display_progress && (mod(step, max(1, floor(n_steps/20))) == 0 || step == n_steps)
        fprintf('    t = %6.0f s (%5.1f min), h_max = %.4f m, Q_out = %.4f m^3/s\n', ...
            t_new, t_new/60, max(d_tot(:))/1000, Q_out_m3s);
    end
end

comp_time = toc;

%% Trim output arrays
t_save = t_save(1:save_idx);
Q_outlet_save = Q_outlet_save(1:save_idx);
if n_track > 0
    h_track_save = h_track_save(:, 1:save_idx);
end

%% Build results struct
results.t = t_save;
results.Q_outlet = Q_outlet_save;
results.h_final = d_tot / 1000;  % Convert mm -> m
results.solver_name = 'CA_Routing';
results.computation_time = comp_time;
if n_track > 0
    results.h_track = h_track_save;       % [n_track x n_save] depth [m] at TrackCells
    results.track_cells = track_cells;    % [n_track x 2] (row, col)
end

if use_infiltration
    results.F_final = F_cum;
else
    results.F_final = zeros(ny, nx);
end

% Water balance
% ORIGINAL (HydroPol2D):
%   total_ponding = sum(results.h_final(:)) * cell_area;
% This summed ponding depth over ALL cells including invalid ones.
if isfield(sim_params, 'valid_mask')
    total_ponding = sum(results.h_final(sim_params.valid_mask)) * cell_area;
else
    total_ponding = sum(results.h_final(:)) * cell_area;
end
results.water_balance.rainfall = total_rainfall_vol;
results.water_balance.infiltration = total_infiltration_vol;
results.water_balance.ponding = total_ponding;
results.water_balance.runoff = total_outflow_vol;
results.water_balance.balance_error = total_rainfall_vol ...
    - total_infiltration_vol - total_ponding - total_outflow_vol;

if display_progress
    fprintf('\n  CA_Routing Results:\n');
    fprintf('    Computation time: %.2f s\n', comp_time);
    fprintf('    Peak outlet Q: %.4f m^3/s\n', max(Q_outlet_save));
    fprintf('    Final h: mean=%.4f m, max=%.4f m\n', ...
        mean(results.h_final(:)), max(results.h_final(:)));
    fprintf('    Water balance:\n');
    fprintf('      Rainfall:     %10.1f m^3\n', total_rainfall_vol);
    fprintf('      Infiltration: %10.1f m^3\n', total_infiltration_vol);
    fprintf('      Ponding:      %10.1f m^3\n', total_ponding);
    fprintf('      Runoff:       %10.1f m^3\n', total_outflow_vol);
    fprintf('      Error:        %10.1f m^3 (%.2f%%)\n', ...
        results.water_balance.balance_error, ...
        100 * abs(results.water_balance.balance_error) / max(total_rainfall_vol, 1e-10));
end

end
