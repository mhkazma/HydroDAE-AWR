%% params_vtilted_AWR2026.m
% Parameter definitions for V-tilted synthetic watershed
% Used in AWR2026 paper for DAE validation and state estimation

% Domain characteristics:
%   - V-shaped cross-section (analytical geometry)
%   - Uniform rainfall input
%   - No infiltration (simplified hydrology)
%   - Synthetic discharge observations
%
% Author: M. H. Kazma

%% Domain Configuration
domain.name = 'vtilted';
domain.type = 'synthetic';

% Physical domain size (fixed - do not change)
% Lx = 1620 chosen so that dx=20 gives nx=81 (odd) naturally
domain.Lx = 1620;            % Domain length in x-direction (meters)
domain.Ly = 1000;            % Domain length in y-direction (meters)

% ---- MESH RESOLUTION SELECTION ----
% Choose cell size: smaller = finer mesh, smoother hydrograph, slower run
% NOTE: Lateral cell count (nx) is forced ODD so the channel center
%       falls exactly on a cell, ensuring perfect V-symmetry.
if ~exist('VTILTED_DX', 'var')
    VTILTED_DX = 20;  % Default 20m cells (set before calling this script to override)
end
domain.dx = VTILTED_DX;
domain.dy = VTILTED_DX;  % Square cells

% Compute grid dimensions from physical size
domain.nx = round(domain.Lx / domain.dx);
domain.ny = round(domain.Ly / domain.dy);

% Force odd lateral cell count for channel symmetry
if mod(domain.nx, 2) == 0
    domain.nx = domain.nx + 1;  % e.g. 80 -> 81
    fprintf('  NOTE: nx adjusted to %d (odd) for channel symmetry\n', domain.nx);
end

% Domain extent (may differ slightly from Lx/Ly due to odd adjustment)
domain.length_x = domain.nx * domain.dx;
domain.length_y = domain.ny * domain.dy;

% V-shaped geometry parameters
domain.v_angle = 45;         % V-angle in degrees
domain.base_slope = 0.02;    % Longitudinal slope (2%) - matches setup_vtilted_domain.m
domain.valley_depth = 5;     % Maximum depth of V-shape (m)

%% Hydraulic Parameters
hydraulic.manning_n = 0.03;  % Manning's roughness coefficient (base/average)
hydraulic.manning_n_hillslope = 0.015;  % Manning's n for hillslopes
hydraulic.manning_n_channel = 0.15;     % Manning's n for channel cells
hydraulic.lateral_slope = 0.05;         % Lateral slope Sx (5%)
hydraulic.channel_width_cells = 3;      % Channel width in cells (center ± 1)
hydraulic.g = 9.81;          % Gravitational acceleration (m/s²)

% Routing method selection
% Options: 'manning_ca' (cellular automata DAE) or 'local_inertial' (momentum-based) or 'semi_implicit' (explicit comparison)
hydraulic.routing_method = 'manning_ca';  % Use Manning CA with BDF Newton for comparison

% Local Inertial formulation type (only used if routing_method = 'local_inertial')
% Options: 'ode' (explicit momentum, all differential) or 'implicit' (implicit momentum, DAE system)
hydraulic.local_inertial_type = 'ode';
% Note: Bounded friction coefficient prevents stiffness at shallow depths

% Semi-implicit solver settings (only used if routing_method = 'semi_implicit')
hydraulic.semi_implicit_dt = 0.5;  % Time step (seconds) for semi-implicit solver

% Numerical stability
hydraulic.min_depth = 1e-6;  % Minimum water depth (m)
hydraulic.max_velocity = 10; % Maximum velocity for stability (m/s)

%% Initial Conditions
initial.depth = 0.0;         % Initial water depth (m) - dry bed (consistent with q=0)
initial.velocity_x = 0.0;    % Initial velocity x-component (m/s)
initial.velocity_y = 0.0;    % Initial velocity y-component (m/s)
initial.infiltration = 0.0;  % Initial infiltration depth (m)

%% Boundary Conditions
% V-tilted uses natural outflow at domain boundaries
bc.type = 'open';            % Open boundary (zero-gradient)
bc.outflow_edge = 'bottom';  % Outflow at bottom edge (y=0)

% Upstream inflow (if needed)
bc.inflow = 0.0;             % No upstream inflow (m³/s)

%% Rainfall Forcing
% Rainfall configuration - choose type below
rainfall = struct();

% --- Option A: Uniform constant rainfall (entire simulation) ---
% rainfall.type = 'uniform_constant';
% rainfall.intensity = 50.0;    % Rainfall intensity (mm/hr)
% rainfall.duration = 1.5;      % Rainfall duration (hours)
% rainfall.start_time = 0.0;    % Start time (hours)

% --- Option B: Rain schedule (to test R=0 periods) ---
rainfall.type = 'uniform_schedule';
rainfall.intensity = 50.0;    % Peak rainfall intensity (mm/hr) - used by plots
rainfall.total_duration = 1.5; % Total simulation window (hours)
rainfall.schedule = [
    0.00 0.10 0.0   % Rain for first 45 minutes at 50 mm/hr
    0.10 0.15 30.0   % Rain for first 45 minutes at 50 mm/hr
    0.15 0.30 50.0 % Rain restarts after 1 hour for 30 minutes at 50 mm/hr
    0.30 0.45 45.0  
    0.45 0.60 0.0   % Rain stops for 15 minutes (R=0)
    0.60 0.90 40.0 % Rain restarts for last 30 minutes at 50 mm/hr
    0.90 1.20 20.0   % Rain continues for 30 minutes at 20 mm/hr
    1.20 1.50 0.0   
];
% Rain stops from 45-60 min, then restarts
% This tests variance propagation with rain-off-rain cycle

% Convert to m/s for simulation (only for uniform_constant)
if strcmp(rainfall.type, 'uniform_constant') && isfield(rainfall, 'intensity')
    rainfall.intensity_ms = rainfall.intensity / 1000 / 3600;  % mm/hr -> m/s
end

%% Infiltration Parameters
% V-tilted: Use Green-Ampt infiltration for hydrologic component
infiltration.model = 'green_ampt';  % Options: 'none', 'green_ampt', 'philip'
infiltration.enabled = true;        % Include infiltration state (F_cum) in model

% Green-Ampt parameters (Higher permeability for testing PSE on pervious soil)
% Sandy soil to see significant infiltration uncertainty
infiltration.Ks_mmh = 15;        % Saturated hydraulic conductivity (mm/hr) - sandy soil
infiltration.Ks = infiltration.Ks_mmh / 3600 / 1000;  % Convert to m/s
infiltration.psi = 0.088;        % Suction head (m) ~88 mm
infiltration.theta_i = 0.15;     % Initial moisture content (-) wetter conditions
infiltration.theta_s = 0.30;     % Saturated moisture content (-) sandy loam
infiltration.delta_theta = infiltration.theta_s - infiltration.theta_i;  % 0.15

%% Time Configuration
time.start = 0.0;             % Start time (hours)
time.end = 1.5;               % End time (hours) - 90 minutes full simulation
time.dt_output = 0.05;        % Output timestep (hours) - 3 minutes

% Convert to seconds for solver
time.start_s = time.start * 3600;
time.end_s = time.end * 3600;
time.dt_output_s = time.dt_output * 3600;

% Solver tolerances (relaxed for stability)
time.RelTol = 1e-5;           % Relative tolerance for ode15i (relaxed)
time.AbsTol = 1e-6;           % Absolute tolerance for ode15i (relaxed)

%% DAE Configuration
% For Manning CA: 5-state DAE [h, qx, qy, vx, vy] or 6-state [h, qx, qy, vx, vy, F]

% Define state components as parameters
dae.n_base_states = 3;           % h, qx, qy (common to all formulations)
dae.n_algebraic_states = 2;      % vx, vy (only for Manning CA)
dae.n_infiltration_states = 1;   % F_cum (if infiltration is used)

if strcmp(hydraulic.routing_method, 'local_inertial')
    % Local Inertial formulation
    if strcmp(hydraulic.local_inertial_type, 'ode')
        dae.formulation = '4state_ode';  % Local Inertial ODE (all differential)
    else
        dae.formulation = '4state_implicit';  % Local Inertial DAE (implicit momentum)
    end
    % Number of states: base (h, qx, qy) + infiltration if enabled
    dae.n_states_per_cell = dae.n_base_states + double(infiltration.enabled) * dae.n_infiltration_states;
else
    % Manning CA formulation
    if infiltration.enabled
        dae.formulation = '6state';  % Manning CA with infiltration
    else
        dae.formulation = '5state';  % Manning CA without infiltration
    end
    % Number of states: base (h, qx, qy) + algebraic (vx, vy) + infiltration if enabled
    dae.n_states_per_cell = dae.n_base_states + dae.n_algebraic_states + double(infiltration.enabled) * dae.n_infiltration_states;
end
dae.solver = 'ode15i';               % MATLAB implicit DAE solver
dae.use_consistent_ic = true;        % Compute consistent initial conditions

% Solver method for Manning CA formulation
% Options: 'bdf' (BDF fixed-point), 'bdf_newton' 
%          'ode15s' (semi-explicit DAE)
dae.solver_method = 'bdf';

% Regularization for ode15s DAE solver (only used when solver_method = 'ode15s')
% Small epsilon converts algebraic constraints to fast dynamics: eps*dQ/dt = Q_manning - Q
dae.epsilon = 0.01;  % 10ms time constant for Q dynamics

% ode15s settings (tolerances for adaptive time-stepping)
dae.RelTol = 1e-4;          % Relative tolerance for ode15s
dae.AbsTol = 1e-6;          % Absolute tolerance for ode15s
dae.MaxStep = 10.0;         % Maximum time step [s]
dae.OutputInterval = 30;    % Output interval [s] (controls tspan density)

% BDF fixed-step settings (used by 'bdf' and 'bdf_newton' solver_methods)
dae.bdf_dt = 1;             % Time step [s] for fixed-step BDF solvers
dae.bdf_order = 1;          % 1 = Backward Euler, 2 = BDF2
dae.bdf_newton_tol = 1e-10;  % Newton iteration tolerance
dae.use_chordal = true;     % Chordal decomposition for bdf_newton solvers

% DAE system size
dae.n_cells = domain.nx * domain.ny;
dae.n_total = dae.n_states_per_cell * dae.n_cells;

%% Observation Network (for State Estimation)
% Synthetic observations at selected cells
obs.n_stations = 5;           % Number of observation points
obs.type = 'discharge';       % Observe discharge (could be depth, velocity)

% Observation locations - defined as fractions of domain, converted to cell indices
% Format: [frac_x, frac_y] where (0,0) = top-left, (1,1) = bottom-right
obs_frac = [
    0.50, 0.20;   % Upstream center
    0.30, 0.50;   % Mid-left
    0.70, 0.50;   % Mid-right
    0.50, 0.80;   % Downstream center
    0.50, 0.96    % Near outlet
];
obs.locations = round([obs_frac(:,1) * domain.nx, obs_frac(:,2) * domain.ny]);
obs.locations = max(obs.locations, 1);  % Ensure >= 1
obs.locations(:,1) = min(obs.locations(:,1), domain.nx);
obs.locations(:,2) = min(obs.locations(:,2), domain.ny);

% Observation uncertainty
obs.std_depth = 0.01;         % Standard deviation for depth (m)
obs.std_discharge = 0.05;     % Standard deviation for discharge (m³/s)
obs.frequency = 0.1;          % Observation frequency (hours) - 6 minutes

%% State Estimation Configuration
estimation.method = 'closed_form';  %  closed-form algorithm
estimation.use_covariance = true;   % Use process covariance matrix

% Process noise parameters
estimation.sigma_h = 0.001;   % Depth process noise std (m)
estimation.sigma_q = 0.01;    % Discharge process noise std (m³/s)
estimation.sigma_v = 0.01;    % Velocity process noise std (m/s)

% Measurement noise (matches observation uncertainty)
estimation.R_depth = obs.std_depth^2;
estimation.R_discharge = obs.std_discharge^2;

%% PSE Distribution Configuration
% Distribution type for state variables (for CI computation)
% Options: 'gaussian' (default), 'lognormal', 'truncated_gaussian'
pse.distribution_type = 'lognormal';  % Use log-normal for positive states

% Distribution-specific parameters
pse.truncation_lower = 0;      % Lower bound for truncated Gaussian (h >= 0)
pse.min_depth_for_gaussian = 1e-3;  % Below this, use transformed distribution (m)

% Confidence level for intervals
pse.confidence_level = 0.95;   % 95% confidence intervals

% Parameter distributions (for uncertainty inputs)
% Defaults: lognormal for positive parameters, truncated Gaussian for bounded
pse.param_distributions.rainfall = 'lognormal';
pse.param_distributions.manning_n = 'lognormal';
pse.param_distributions.Ks = 'lognormal';
pse.param_distributions.psi = 'lognormal';
pse.param_distributions.eta = 'truncated_gaussian';

% Parameter bounds used when truncated Gaussian is selected
pse.param_bounds.rainfall = [0, inf];
pse.param_bounds.manning_n = [0, inf];
pse.param_bounds.Ks = [0, inf];
pse.param_bounds.psi = [0, inf];
pse.param_bounds.eta = [0, 1];

%% Output Configuration
output.directory = 'results/vtilted/AWR2026';
output.save_states = true;     % Save full state variables
output.save_diagnostics = true; % Save solver diagnostics
output.save_observations = true; % Save synthetic observations
output.save_estimates = true;   % Save state estimates

% Variables to save
output.variables = {'depth', 'discharge_x', 'discharge_y', ...
                    'velocity_x', 'velocity_y'};

% Spatial outputs
output.save_spatial_maps = true;  % Save 2D spatial maps
output.map_times = [0.5, 1.0, 1.5, 2.0];  % Times to save maps (hours)

%% Visualization Settings
viz.create_plots = true;
viz.plot_study_area = true;  % Study area visualization
viz.plot_pse_uncertainty_bounds = true;  % PSE uncertainty bounds (time series)
viz.plot_pse_spatial_uncertainty = true; % PSE spatial uncertainty maps
viz.plot_hydrograph = true;    % Hydrograph at outlet
viz.plot_spatial = true;       % Spatial maps at key times
viz.plot_comparison = true;    % DAE vs HydroPol2D
viz.plot_estimation = true;    % State estimation results

% Figure export settings
viz.export_format = 'png';     % Options: 'png', 'pdf', 'eps'
viz.export_dpi = 300;          % Resolution for raster formats
viz.export_to_latex = false;   % Export using LaTeX utilities

%% Validation Settings
validation.run_hydropol2d = false;  % Run original HydroPol2D for comparison
validation.compare_solvers = false;  % Compare solvers
validation.compute_metrics = true;   % Compute error metrics (RMSE, NSE, etc.)

% Metrics to compute
validation.metrics = {'RMSE', 'NSE', 'MAE', 'PeakError', 'VolumeError'};

%% Display Configuration Summary
fprintf('\n--- V-Tilted Domain Configuration ---\n');
fprintf('Domain size: %d x %d cells (%.0f x %.0f m)\n', ...
    domain.nx, domain.ny, domain.length_x, domain.length_y);
fprintf('Cell resolution: %.1f x %.1f m\n', domain.dx, domain.dy);
fprintf('V-angle: %.0f degrees, Base slope: %.3f\n', ...
    domain.v_angle, domain.base_slope);
fprintf('Manning n: %.3f\n', hydraulic.manning_n);
fprintf('Infiltration: %s (Ks=%.2e m/s)\n', ...
    infiltration.model, infiltration.Ks);
if isfield(rainfall, 'intensity') && isfield(rainfall, 'duration')
    fprintf('Rainfall: %.0f mm/hr for %.1f hours\n', ...
        rainfall.intensity, rainfall.duration);
elseif isfield(rainfall, 'intensity') && isfield(rainfall, 'total_duration')
    fprintf('Rainfall: %s, %.0f mm/hr peak, %.1f hours total\n', ...
        rainfall.type, rainfall.intensity, rainfall.total_duration);
elseif isfield(rainfall, 'schedule')
    fprintf('Rainfall: schedule with %d pulses\n', size(rainfall.schedule, 1));
end
fprintf('Simulation time: %.1f to %.1f hours\n', time.start, time.end);
fprintf('DAE formulation: %s (%d states per cell)\n', ...
    dae.formulation, dae.n_states_per_cell);
fprintf('Total DAE system size: %d states\n', dae.n_total);
fprintf('Observation stations: %d\n', obs.n_stations);
fprintf('------------------------------------\n');
