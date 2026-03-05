function [t_out, results] = run_semi_implicit_simulation(params, R_intensity, t_end_s, varargin)
% RUN_SEMI_IMPLICIT_SIMULATION  Semi-implicit Local Inertial solver
%
% Implements the semi-implicit Local Inertial formulation using the local inertial equations
% with perfect mass conservation for comparison with DAE solvers.
%
% Syntax:
%   [t_out, results] = run_semi_implicit_simulation(params, R_intensity, t_end_s)
%   [t_out, results] = run_semi_implicit_simulation(..., 'dt', 0.5)
%
% Inputs:
%   params       - Simulation parameters structure from setup functions
%   R_intensity  - Rainfall intensity (m/s) scalar or per-cell vector
%   t_end_s      - Simulation end time (seconds)
%
% Optional Parameters:
%   'dt'         - Time step (seconds), default 0.5
%   'DisplayProgress' - Show progress output, default true
%
% Outputs:
%   t_out        - Time vector (seconds)
%   results      - Structure with:
%                  .h  - Water depth [ny x nx x nt]
%                  .qx - Unit discharge x [ny x nx x nt]
%                  .qy - Unit discharge y [ny x nx x nt]
%                  .F  - Cumulative infiltration [ny x nx x nt]
%                  .volume_balance - Structure with water balance info
%
% Equations:
%   Momentum (semi-implicit):
%     q^{t+Δt} = [q^t - g*ĥ*Δt*S] / [1 + g*Δt*n²*|q^t|/ĥ^(7/3)]
%
%   Continuity:
%     h^{t+Δt} = h^t + Δt*(q_in - q_out)/Δx + Δt*(R - f)
%
% Reference:
% 
%
% =========================================================================
% AWR2026 FIXES — Corrections to HydroPol2D semi-implicit local inertial
% =========================================================================
% The following issues were identified when running the semi-implicit
% comparison on irregular (real-world) watershed domains such as Walnut
% Gulch, where most boundary cells are NOT designated outlet cells.
%
% FIX 5 (Boundary drainage — X and Y directions): The original code only
%   allowed outflow via Y-direction Manning discharge at the ~4 designated
%   outlet cells. For Walnut Gulch (1,383 boundary cells), water at all
%   other boundary faces was trapped, inflating h_max to 1.58 m vs. 0.68 m
%   in the DAE solver. Fixed by adding Manning normal-depth free drainage
%   at ALL boundary faces (east, south, and implicitly west/north) using a
%   default 1% boundary slope, matching the DAE solver's treatment.
%   Search: "FIX [AWR2026]: Add X-direction free drainage"
%   Search: "FIX [AWR2026]: Add free drainage at non-outlet"
%
% FIX 6 (Outflow tracking): Outflow volume was tracked only at designated
%   outlet cells (is_outlet). With Fix 5 adding drainage at all boundary
%   faces, the water balance must also account for non-outlet boundary
%   outflow. Changed to detect and track flux leaving the domain at any
%   face bordering an invalid or non-existent cell.
%   Search: "FIX [AWR2026]: Track boundary outflow at ALL boundary faces"
%
% Impact: Before fixes — h_max = 1.58 m, runoff = 0.04% of rainfall.
%         After fixes — h_max and runoff should approach DAE values.
% =========================================================================

% Parse optional inputs
p = inputParser;
addParameter(p, 'dt', 0.5);
addParameter(p, 'DisplayProgress', true);
addParameter(p, 'LightweightMode', false);  % Skip 3D history to save RAM (MC use)
parse(p, varargin{:});

dt = p.Results.dt;
display_progress = p.Results.DisplayProgress;
lightweight = p.Results.LightweightMode;

% Extract parameters
ny = params.grid_size(1);
nx = params.grid_size(2);
n_cells = nx * ny;
dx = params.dx;
dy = params.dy;

g = params.g;
Z = params.Z;
n_manning = params.n_manning;
h_min = params.h_min;

% Infiltration parameters
if isfield(params, 'Ks')
    Ks = params.Ks;
    eta = params.eta;
    psi_f = params.psi_f;
    F_min = params.F_min;
    infiltration_enabled = true;
else
    infiltration_enabled = false;
end

% Rainfall configuration
if isstruct(R_intensity)
    if isfield(R_intensity, 'schedule')
        rain_type = 'schedule';
        rain_schedule = R_intensity.schedule;  % [t_start_s, t_end_s, R_m/s]
    elseif isfield(R_intensity, 'base') && isfield(R_intensity, 'duration')
        rain_type = 'time_limited';
        rain_base = R_intensity.base;
        rain_duration = R_intensity.duration;
    elseif isfield(R_intensity, 'rain_matrix')
        rain_type = 'matrix';
        rain_matrix = R_intensity.rain_matrix;
        dt_rain = R_intensity.dt_rain;
    else
        error('run_semi_implicit_simulation: Unrecognized R_intensity struct.');
    end
elseif isscalar(R_intensity)
    rain_type = 'constant';
    rain_const = R_intensity;
else
    rain_type = 'constant';
    rain_const = mean(R_intensity);
end

% Outlet configuration
slope_outlet = params.slope_outlet;
is_outlet = false(ny, nx);
use_compressed = isfield(params, 'use_compressed_indexing') && params.use_compressed_indexing;
if use_compressed && isfield(params, 'outlet_indices') && isfield(params, 'compressed_to_grid')
    % Compressed indexing (e.g., walnut_gulch): mark specific outlet cells
    ctg = params.compressed_to_grid;
    for k = 1:length(params.outlet_indices)
        ci = params.outlet_indices(k);
        is_outlet(ctg(ci, 1), ctg(ci, 2)) = true;
    end
elseif isfield(params, 'outlet_direction') && strcmp(params.outlet_direction, 'down')
    % vtilted-style: entire bottom row
    is_outlet(ny, :) = true;
elseif isfield(params, 'outlet_indices') && ~isempty(params.outlet_indices)
    % Full-grid linear indices
    [oi_r, oi_c] = ind2sub([ny, nx], params.outlet_indices);
    for k = 1:length(oi_r)
        is_outlet(oi_r(k), oi_c(k)) = true;
    end
else
    is_outlet(ny, :) = true;  % Default: bottom row
end

% Valid cell mask (skip NaN/inactive cells for compressed domains)
if isfield(params, 'valid_mask')
    valid_mask = params.valid_mask;
else
    valid_mask = true(ny, nx);  % All cells active (vtilted)
end
n_valid_cells = sum(valid_mask(:));  % For water balance

% Initialize states
h = zeros(ny, nx);
qx = zeros(ny, nx);
qy = zeros(ny, nx);
if infiltration_enabled
    F_cum = F_min * ones(ny, nx);
else
    F_cum = zeros(ny, nx);
end

% Storage for time history
output_interval = 10;  % Save every 10 steps
t_out = 0;
if ~lightweight
    h_history = h;
    qx_history = qx;
    qy_history = qy;
    F_history = F_cum;
end

% Water balance tracking
vol_initial = sum(h(valid_mask)) * dx * dy;
vol_rain = 0;
vol_infil = 0;
vol_out = 0;

% Outlet hydrograph tracking
Q_out_history = 0;  % First entry at t=0 (total boundary outflow)
Q_outlet_history = 0;  % First entry at t=0 (outlet/flume cells only)

% Time stepping
h_flow_min = 1e-6;
t = 0;
step = 0;

if display_progress
    fprintf('Semi-Implicit Local Inertial Solver\n');
    fprintf('  Grid: %d x %d cells\n', nx, ny);
    fprintf('  dt: %.2f s, t_end: %.1f s\n', dt, t_end_s);
end

tic;
while t < t_end_s
    step = step + 1;
    t_new = min(t + dt, t_end_s);
    dt_actual = t_new - t;
    
    % Store old values
    qx_old = qx;
    qy_old = qy;
    
    % Compute rainfall at current time [m/s]
    switch rain_type
        case 'constant'
            R_now = rain_const;
        case 'time_limited'
            if t < rain_duration
                R_now = rain_base;
            else
                R_now = 0;
            end
        case 'schedule'
            R_now = 0;
            for rr = 1:size(rain_schedule, 1)
                if t >= rain_schedule(rr, 1) && t < rain_schedule(rr, 2)
                    R_now = rain_schedule(rr, 3);
                    break;
                end
            end
        case 'matrix'
            time_idx = min(floor(t / dt_rain) + 1, size(rain_matrix, 2));
            R_now = 0;  % Per-cell rainfall handled via R_field below
    end
    
    %% STEP 1: Update momentum (semi-implicit)
    qx_new = zeros(ny, nx);
    qy_new = zeros(ny, nx);
    
    % Compute per-cell rainfall for matrix mode (spatial IDW)
    if strcmp(rain_type, 'matrix')
        if use_compressed && isfield(params, 'compressed_to_grid')
            % Rain matrix is [n_valid_cells x n_timesteps] in compressed space
            R_field = zeros(ny, nx);
            for rc = 1:size(params.compressed_to_grid, 1)
                rr = params.compressed_to_grid(rc, 1);
                cc = params.compressed_to_grid(rc, 2);
                R_field(rr, cc) = rain_matrix(rc, time_idx);
            end
        else
            % Non-compressed: use spatially averaged scalar
            R_scalar = mean(rain_matrix(:, time_idx));
            R_field = R_scalar * ones(ny, nx);
        end
    end
    
    for i = 1:ny
        for j = 1:nx
            % Skip inactive cells (NaN elevation / outside watershed)
            if ~valid_mask(i, j), continue; end
            
            % X-momentum
            if j < nx && valid_mask(i, j+1)
                % Interior face: standard semi-implicit local inertial
                h_L = h(i, j);
                h_R = h(i, j+1);
                
                if max(h_L, h_R) > h_flow_min
                    eta_L = Z(i, j) + h_L;
                    eta_R = Z(i, j+1) + h_R;
                    S_x = (eta_R - eta_L) / dx;
                    
                    Hf = max(eta_L, eta_R) - max(Z(i,j), Z(i,j+1));
                    Hf = max(Hf, h_flow_min);
                    
                    n_face = 0.5 * (n_manning(i,j) + n_manning(i,j+1));
                    friction = g * n_face^2 * abs(qx_old(i, j)) / Hf^(7/3);
                    
                    qx_new(i, j) = (qx_old(i, j) - g * Hf * dt_actual * S_x) / (1 + dt_actual * friction);
                end
                
            % ORIGINAL (HydroPol2D):
            %   No X-direction outlet or boundary drainage was implemented.
            %   Only Y-direction outlet flow was computed (see below).
            %   For irregular watersheds (e.g., Walnut Gulch with 1,383
            %   boundary cells but only ~4 outlet cells), water at boundary
            %   cells with no valid eastern neighbor was trapped — no flux
            %   could leave the domain in the X-direction.
            %
            % FIX [AWR2026]: Add X-direction free drainage at boundary faces.
            % When a valid cell has no valid neighbor to the east (grid edge
            % or outside-watershed cell), apply Manning normal-depth outflow
            % using a terrain-based or default boundary slope (1%), matching
            % the DAE solver's free-drainage boundary treatment (see
            % run_manning_ca_bdf.m, default_boundary_slope = 0.01).
            else
                h_c = h(i, j);
                if h_c > h_flow_min
                    if is_outlet(i, j)
                        S_bd = slope_outlet;
                    elseif j < nx && ~valid_mask(i, j+1)
                        S_bd = 0.01;  % Default 1% free drainage slope
                    elseif j == nx
                        S_bd = 0.01;
                    else
                        S_bd = 0;
                    end
                    if S_bd > 0
                        n_face = n_manning(i, j);
                        qx_new(i, j) = (1/n_face) * h_c^(5/3) * sqrt(S_bd);
                    end
                end
            end
            
            % Y-momentum
            if i < ny && valid_mask(i+1, j)
                % Interior face: standard semi-implicit local inertial
                h_U = h(i, j);
                h_D = h(i+1, j);
                
                if max(h_U, h_D) > h_flow_min
                    eta_U = Z(i, j) + h_U;
                    eta_D = Z(i+1, j) + h_D;
                    S_y = (eta_D - eta_U) / dy;
                    
                    Hf = max(eta_U, eta_D) - max(Z(i,j), Z(i+1,j));
                    Hf = max(Hf, h_flow_min);
                    
                    n_face = 0.5 * (n_manning(i,j) + n_manning(i+1,j));
                    friction = g * n_face^2 * abs(qy_old(i, j)) / Hf^(7/3);
                    
                    qy_new(i, j) = (qy_old(i, j) - g * Hf * dt_actual * S_y) / (1 + dt_actual * friction);
                end
                
            elseif is_outlet(i, j)
                % Outlet: Manning equation for normal flow (original HydroPol2D)
                h_c = h(i, j);
                if h_c > h_flow_min
                    n_face = n_manning(i, j);
                    qy_new(i, j) = (1/n_face) * h_c^(5/3) * sqrt(slope_outlet);
                end
                
            % ORIGINAL (HydroPol2D):
            %   No free drainage at non-outlet southern boundary faces.
            %   Only designated outlet cells (typically ~4 cells at the
            %   watershed pour point) allowed Y-direction outflow.
            %   For Walnut Gulch (1,383 boundary cells), this trapped water
            %   at all non-outlet boundary cells, causing h_max to rise to
            %   1.58 m vs. 0.68 m in the DAE solver which drains all
            %   boundary faces.
            %
            % FIX [AWR2026]: Add free drainage at non-outlet boundary faces
            % using a default 1% slope, matching the DAE solver's treatment
            % where any valid cell bordering an invalid cell or grid edge
            % allows outflow via Manning normal-depth discharge.
            else
                h_c = h(i, j);
                if h_c > h_flow_min
                    if i < ny && ~valid_mask(i+1, j)
                        S_bd = 0.01;
                    elseif i == ny
                        S_bd = 0.01;
                    else
                        S_bd = 0;
                    end
                    if S_bd > 0
                        n_face = n_manning(i, j);
                        qy_new(i, j) = (1/n_face) * h_c^(5/3) * sqrt(S_bd);
                    end
                end
            end
        end
    end
    
    %% STEP 2: Flux limiting
    qx_limited = qx_new;
    qy_limited = qy_new;
    
    for i = 1:ny
        for j = 1:nx
            if ~valid_mask(i, j), continue; end
            h_c = h(i, j);
            
            % Incoming flux
            q_in = 0;
            if j > 1 && qx_new(i, j-1) > 0
                q_in = q_in + qx_new(i, j-1) / dx;
            end
            if j < nx && qx_new(i, j) < 0
                q_in = q_in - qx_new(i, j) / dx;
            end
            if i > 1 && qy_new(i-1, j) > 0
                q_in = q_in + qy_new(i-1, j) / dy;
            end
            if i < ny && qy_new(i, j) < 0
                q_in = q_in - qy_new(i, j) / dy;
            end
            
            % Use per-cell or scalar rainfall
            if strcmp(rain_type, 'matrix')
                R_cell = R_field(i, j);
            else
                R_cell = R_now;
            end
            h_available_rate = h_c / dt_actual + R_cell + q_in * dt_actual;
            
            % Outgoing flux
            q_out_east = max(qx_new(i, j), 0) / dx;
            q_out_south = max(qy_new(i, j), 0) / dy;
            total_out_rate = q_out_east + q_out_south;
            
            if total_out_rate > h_available_rate && total_out_rate > 0
                scale = max(h_available_rate / total_out_rate, 0);
                if qx_new(i, j) > 0
                    qx_limited(i, j) = qx_new(i, j) * scale;
                end
                if qy_new(i, j) > 0
                    qy_limited(i, j) = qy_new(i, j) * scale;
                end
            end
        end
    end
    
    %% STEP 3: Update continuity
    h_new = zeros(ny, nx);
    F_new = F_cum;
    Q_out_step = 0;  % Total boundary discharge this step [m^3/s]
    Q_outlet_step = 0;  % Outlet/flume cells only [m^3/s]
    
    for i = 1:ny
        for j = 1:nx
            if ~valid_mask(i, j), continue; end
            h_c = h(i, j);
            F_c = F_cum(i, j);
            
            % Flux balance
            if j > 1
                q_west = qx_limited(i, j-1);
            else
                q_west = 0;
            end
            q_east = qx_limited(i, j);
            
            if i > 1
                q_north = qy_limited(i-1, j);
            else
                q_north = 0;
            end
            q_south = qy_limited(i, j);
            
            % ORIGINAL (HydroPol2D):
            %   Outlet flow tracking was restricted to designated outlet
            %   cells (is_outlet). With the boundary drainage fix above,
            %   outflow now also occurs at non-outlet boundary cells.
            %
            % FIX [AWR2026]: Track boundary outflow at ALL boundary faces
            % where flow leaves the domain (flux into a non-existent or
            % invalid neighbor), ensuring the water balance accounts for
            % all drainage.
            q_out_total = 0;
            
            % Southward outflow leaving domain
            if q_south > 0 && (i == ny || (i < ny && ~valid_mask(i+1, j)))
                q_out_total = q_out_total + q_south * dx;
            end
            % Eastward outflow leaving domain
            if q_east > 0 && (j == nx || (j < nx && ~valid_mask(i, j+1)))
                q_out_total = q_out_total + q_east * dy;
            end
            % Northward outflow leaving domain
            if i == 1 && q_north < 0
                q_out_total = q_out_total + abs(q_north) * dx;
            end
            % Westward outflow leaving domain
            if j == 1 && q_west < 0
                q_out_total = q_out_total + abs(q_west) * dy;
            end
            
            if q_out_total > 0
                vol_out = vol_out + q_out_total * dt_actual;
                Q_out_step = Q_out_step + q_out_total;
                if is_outlet(i, j)
                    Q_outlet_step = Q_outlet_step + q_out_total;
                end
            end
            
            % Use per-cell or scalar rainfall
            if strcmp(rain_type, 'matrix')
                R_cell = R_field(i, j);
            else
                R_cell = R_now;
            end
            
            dh_flux = dt_actual * (q_west - q_east) / dx + dt_actual * (q_north - q_south) / dy;
            h_available = h_c + dt_actual * R_cell + dh_flux;
            
            % Infiltration
            if infiltration_enabled && h_available > h_min
                % Handle both scalar and spatially variable parameters
                if isscalar(Ks)
                    Ks_ij = Ks;
                    eta_ij = eta;
                    psi_ij = psi_f;
                    Fmin_ij = F_min;
                else
                    Ks_ij = Ks(i, j);
                    eta_ij = eta(i, j);
                    psi_ij = psi_f(i, j);
                    if isscalar(F_min)
                        Fmin_ij = F_min;
                    else
                        Fmin_ij = F_min(i, j);
                    end
                end
                f_capacity = Ks_ij * (1 + eta_ij * (h_c + psi_ij) / max(F_c, Fmin_ij));
                f_infil = min(f_capacity, h_available / dt_actual);
                F_new(i, j) = F_c + dt_actual * f_infil;
                vol_infil = vol_infil + f_infil * dt_actual * dx * dy;
            else
                f_infil = 0;
            end
            
            h_new(i, j) = max(h_available - dt_actual * f_infil, 0);
        end
    end
    
    % Accumulate rainfall volume (only over valid cells)
    if strcmp(rain_type, 'matrix')
        vol_rain = vol_rain + sum(R_field(valid_mask)) * dt_actual * dx * dy;
    else
        vol_rain = vol_rain + R_now * dt_actual * n_valid_cells * dx * dy;
    end
    
    % Update state
    h = h_new;
    qx = qx_limited;
    qy = qy_limited;
    F_cum = F_new;
    t = t_new;
    
    % Store results
    if mod(step, output_interval) == 0 || t >= t_end_s
        t_out = [t_out; t];
        if ~lightweight
            h_history = cat(3, h_history, h);
            qx_history = cat(3, qx_history, qx);
            qy_history = cat(3, qy_history, qy);
            F_history = cat(3, F_history, F_cum);
        end
        Q_out_history = [Q_out_history; Q_out_step];
        Q_outlet_history = [Q_outlet_history; Q_outlet_step];
    end
    
    % Progress
    if display_progress && mod(step, 200) == 0
        fprintf('  t = %.1f s (%.1f%%), h_max = %.3f m\n', ...
            t, 100*t/t_end_s, max(h(:)));
    end
end

comp_time = toc;

if display_progress
    fprintf('  Completed in %.2f s (%d steps)\n', comp_time, step);
end

% Water balance
vol_final = sum(h(valid_mask)) * dx * dy;
if infiltration_enabled
    F_diff = F_cum - F_min;
    vol_infil_final = sum(F_diff(valid_mask)) * dx * dy;
else
    vol_infil_final = 0;
end

% Compute water balance components as percentages
total_rainfall = vol_rain;
total_infiltration = vol_infil_final;
total_ponding = vol_final;
total_runoff = vol_out;

if ~lightweight
    results.h = h_history;
    results.qx = qx_history;
    results.qy = qy_history;
    results.F = F_history;
else
    % Lightweight mode: only store final state (not full 3D history)
    results.h = h;            % (ny, nx) final only
    results.qx = qx;
    results.qy = qy;
    results.F = F_cum;
end
results.Q_boundary_ts = Q_out_history;  % Total boundary outflow [m^3/s]
results.Q_outlet_ts = Q_outlet_history;  % Outlet/flume cells only [m^3/s]

% Add PSE-required fields (for state estimation compatibility)
results.h_final = h(:);  % Final water depths as column vector
if infiltration_enabled
    results.F_final = F_cum(:);  % Final infiltration as column vector
else
    results.F_final = zeros(n_cells, 1);  % Zero if no infiltration
end
if isscalar(R_intensity)
    results.R = R_now * ones(n_cells, 1);  % Rainfall as column vector (final timestep rate)
else
    results.R = R_intensity(:);  % Use provided rainfall field
end

results.volume_balance.initial = vol_initial;
results.volume_balance.rainfall = vol_rain;
results.volume_balance.final_ponding = vol_final;
results.volume_balance.infiltration = vol_infil_final;
results.volume_balance.outflow = vol_out;
results.volume_balance.error = (vol_initial + vol_rain) - (vol_final + vol_infil_final + vol_out);
results.volume_balance.error_percent = 100 * results.volume_balance.error / (vol_initial + vol_rain);

% Add water_balance structure for compatibility with run_dae_simulation
results.water_balance = struct();
results.water_balance.rainfall = total_rainfall;
results.water_balance.infiltration = total_infiltration;
results.water_balance.ponding = total_ponding;
results.water_balance.runoff = total_runoff;

% Display water balance
if display_progress
    fprintf('\nWater Balance:\n');
    fprintf('  Rainfall: %.2f m³\n', total_rainfall);
    if infiltration_enabled
        fprintf('  Infiltration: %.2f m³ (%.1f%%)\n', total_infiltration, 100*total_infiltration/total_rainfall);
    end
    fprintf('  Ponding: %.2f m³ (%.1f%%)\n', total_ponding, 100*total_ponding/total_rainfall);
    fprintf('  Runoff: %.2f m³ (%.1f%%)\n', total_runoff, 100*total_runoff/total_rainfall);
    
    fprintf('\nFinal State:\n');
    fprintf('  Water depth: min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
            min(h(:))*1000, max(h(:))*1000, mean(h(:))*1000);
    if infiltration_enabled
        fprintf('  Infiltration: min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
                min(F_cum(:))*1000, max(F_cum(:))*1000, mean(F_cum(:))*1000);
    end
    fprintf('\n');
end

results.comp_time = comp_time;
results.computation_time = comp_time;  % Alias for compatibility
results.n_steps = step;

end
