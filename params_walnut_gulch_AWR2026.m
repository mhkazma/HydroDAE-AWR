%% params_walnut_gulch_AWR2026.m
% Parameter definitions for Walnut Gulch Experimental Watershed
% Used in AWR2026 paper for real-world DAE application and state estimation
%
% Watershed characteristics:
%   - Location: Tombstone, Arizona, USA
%   - Area: ~149 km² (full watershed)
%   - Climate: Semi-arid with monsoon rainfall
%   - Land use: Rangeland, minimal vegetation
%   - Well-instrumented with rain and stream gauges
%
% Author: M. H. Kazma

%% Domain Configuration
domain.name = 'walnut_gulch';
domain.type = 'real_world';

% Grid dimensions - Based on target resolution and watershed extent
% Available grids (created by create_walnut_gulch_grid.m):
domain.grid_resolution = [242, 521];  % [nrows, ncols] - 50m resolution
domain.nx = domain.grid_resolution(2);  % Number of cells in x-direction (cols)
domain.ny = domain.grid_resolution(1);  % Number of cells in y-direction (rows)

% Cell size (computed from watershed extent and grid resolution)
domain.watershed_extent_x = 26286;  % E-W extent (m) from watershed bounds
domain.watershed_extent_y = 12596;  % N-S extent (m) from watershed bounds
domain.dx = domain.watershed_extent_x / domain.nx;  % Cell size x (~100.7 m)
domain.dy = domain.watershed_extent_y / domain.ny;  % Cell size y (~104.1 m)

% Domain extent
domain.area_km2 = 149;       % Watershed area (km²)
domain.outlet_name = 'Flume 1';  % Main outlet gauge

% Processed data files (created by data processing scripts)
domain.spatial_data_file = 'data/walnut_gulch/processed/walnut_gulch_spatial_data.mat';
% Grid file naming: {nrows}x{ncols} format
% Use stream-burned DEM for hydrological connectivity
domain.use_burned_dem = true;  % Set to false to use original DEM
if domain.use_burned_dem
    domain.grid_file = sprintf('data/walnut_gulch/processed/walnut_gulch_grid_%dx%d_burned.mat', ...
                               domain.grid_resolution(1), domain.grid_resolution(2));
else
    domain.grid_file = sprintf('data/walnut_gulch/processed/walnut_gulch_grid_%dx%d.mat', ...
                               domain.grid_resolution(1), domain.grid_resolution(2));
end

%% Hydraulic Parameters
% Routing method selection
% Options: 'manning_ca' (cellular automata DAE) or 'local_inertial' (momentum-based) or 'semi_implicit' (explicit comparison)
hydraulic.routing_method = 'manning_ca';  % Semi-implicit routing method

% Local Inertial formulation type (only used if routing_method = 'local_inertial')
% Options: 'ode' (explicit momentum, all differential) or 'implicit' (implicit momentum, DAE system)
hydraulic.local_inertial_type = 'ode';

% Semi-implicit solver settings (only used if routing_method = 'semi_implicit')
hydraulic.semi_implicit_dt = 0.5;  % Time step (seconds) for semi-implicit solver

% Manning's roughness - spatially variable from 24 vegetation types (IMPLEMENTED)
hydraulic.manning_n_default = 0.035;  % Default for rangeland (fallback only)
hydraulic.manning_n_min = 0.015;      % Minimum (built-up, bare soil)
hydraulic.manning_n_max = 0.078;      % Maximum (trees and shrub along drains)
hydraulic.use_spatial_manning = true; % Grid has spatially variable Manning's n

% Manning's n spatial statistics (from vegetation mapping)
% Mapped from 24 vegetation cover types
% Result: 5.0x spatial variation (0.015 - 0.075)
% Coverage: 100% of grid cells (291/291)
% Data file: data/walnut_gulch/processed/vegetation_parameters.csv

% Physical constants
hydraulic.g = 9.81;          % Gravitational acceleration (m/s²)

% Numerical stability
hydraulic.min_depth = 1e-6;  % Minimum water depth (m)
hydraulic.max_velocity = 15; % Maximum velocity for stability (m/s)

%% Initial Conditions
initial.depth = 0.0;         % Initial water depth (m) - dry bed (pre-storm)
initial.velocity_x = 0.0;    % Initial velocity x-component (m/s)
initial.velocity_y = 0.0;    % Initial velocity y-component (m/s)
initial.infiltration = 0.0;  % Initial infiltration depth (m)
initial.soil_moisture = 0.15; % Initial soil moisture content (-)

%% Boundary Conditions
bc.type = 'open';            % Open boundary at outlet
bc.outflow_location = 'outlet_cells';  % Cells identified as outlet

% Outlet cells (to be identified from flow accumulation)
bc.outlet_cells = [];        % Will be computed from DEM

% Upstream inflow (typically none for closed watershed)
bc.inflow = 0.0;             % No external inflow (m³/s)

%% Rainfall Forcing
% Walnut Gulch supports three rainfall modes:
%   'uniform_constant' - Constant uniform rainfall (for testing / synthetic runs)
%   'uniform_schedule'  - Multi-pulse schedule, spatially uniform
%   'real_event'        - Observed breakpoint data from USDA-ARS gauge network
%
% >>> Change this line to switch rainfall mode <<<
rainfall = struct();
rainfall.type = 'uniform_schedule';  % Options: 'uniform_constant' | 'uniform_schedule' | 'real_event'

% ---- Mode 1: Uniform constant rainfall (synthetic / quick testing) ----
rainfall.intensity = 40;         % Rainfall intensity (mm/hr)
rainfall.duration = 0.5;         % Storm duration (hours) — rain stops, recession follows
rainfall.start_time = 0.0;       % Start time (hours)
% Also store under legacy field names for backward compatibility
rainfall.uniform_intensity = rainfall.intensity;
rainfall.duration_hours = rainfall.duration;

% ---- Mode 2: Multi-pulse uniform schedule (synthetic) ----
% To use, set rainfall.type = 'uniform_schedule'
rainfall.schedule = [
    0.00 0.15 00.0    % start_hr  end_hr  intensity_mmhr
    0.15 0.50 40.0
    0.50 0.75 25.0
    0.75 1.00  0.0    % recession
];

% ---- Mode 3: Real observed event (USDA-ARS gauge network) ----
% Aug 16, 2022 monsoon storm — 125 gages, 661 timesteps at 1-min resolution
% Gages with shapefile coordinates map to 86 grid cells; remaining 39 gages
% (newer installations without coords) contribute to watershed-average only.
rainfall.event_id = '20220823';  % Event identifier
rainfall.event_data_file = 'data/walnut_gulch/processed/event_20220823/rainfall_timeseries_20220823.csv';
rainfall.event_description = 'Aug 23, 2022 monsoon — 126 gages, 32 mm max, 5 flumes, light event';

% Event sub-window (simulate a portion of the event)
% Set event_start_min = 0 and event_duration_min = Inf for the full event.
% Examples:
%   event_start_min=0,   event_duration_min=60   → first hour
%   event_start_min=120, event_duration_min=60   → minutes 120-180
%   event_start_min=0,   event_duration_min=Inf  → full event (default)
rainfall.event_start_min = 0;       % Offset into event CSV [minutes]
rainfall.event_duration_min = Inf;  % Window length [minutes] (Inf = full event)

% Spatial interpolation settings (used for real_event mode)
rainfall.interpolation = 'IDW';      % Inverse Distance Weighting
rainfall.idw_power = 2;              % IDW power parameter
rainfall.use_watershed_avg = false;  % true = use watershed-avg column (Option A)
                                     % false = IDW from individual gages (Option B)

% Rain gage network — 86 gages with grid cell locations
% Mapped from USDA-ARS digital raingages (NAD83 UTM Zone 12N)
% Source: data/walnut_gulch/raw/Files_SWRC/instrumentation/raingage_digital/
rainfall.n_gages = 86;
rainfall.gage_event_nums = [2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, ...
    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, ...
    35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, ...
    53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, ...
    71, 72, 74, 76, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 100, 384, 398];
rainfall.gage_locations = [...  % [row, col] — 1-indexed, row 1=north (50m grid)
      74,   21;  % RG002
     119,   21;  % RG003
      64,   55;  % RG004
      99,   50;  % RG005
      58,   86;  % RG007
      87,   73;  % RG008
     127,   89;  % RG009
     147,   82;  % RG010
      62,  115;  % RG011
      87,  104;  % RG012
     111,  119;  % RG013
     171,  106;  % RG014
      51,  146;  % RG015
      81,  134;  % RG016
     110,  157;  % RG017
     153,  131;  % RG018
     180,  154;  % RG019
     216,  146;  % RG020
      47,  179;  % RG021
      79,  164;  % RG022
      91,  179;  % RG023
     138,  170;  % RG024
     177,  193;  % RG025
     204,  172;  % RG026
      78,  205;  % RG027
     115,  209;  % RG028
     143,  233;  % RG029
     195,  227;  % RG030
      41,  245;  % RG031
      81,  245;  % RG032
     104,  232;  % RG033
     165,  216;  % RG034
     176,  230;  % RG035
     215,  223;  % RG036
     193,  263;  % RG037
      44,  291;  % RG038
      89,  266;  % RG039
     109,  264;  % RG040
     153,  262;  % RG041
     216,  250;  % RG042
      26,  304;  % RG043
      82,  296;  % RG044
     109,  291;  % RG045
     142,  302;  % RG046
     174,  282;  % RG047
     201,  308;  % RG048
     210,  286;  % RG049
      17,  334;  % RG050
      53,  309;  % RG051
      76,  327;  % RG052
     124,  317;  % RG053
      35,  361;  % RG054
      43,  395;  % RG055
      68,  363;  % RG056
      99,  319;  % RG057
     124,  345;  % RG058
     146,  363;  % RG059
      66,  390;  % RG060
      83,  385;  % RG061
      94,  407;  % RG062
     119,  386;  % RG063
      26,  422;  % RG064
      62,  422;  % RG065
      86,  437;  % RG066
      19,  447;  % RG067
      51,  458;  % RG068
       6,  475;  % RG069
      31,  482;  % RG070
      63,  233;  % RG071
      91,  364;  % RG072
      35,  208;  % RG074
     121,   49;  % RG076
     132,  127;  % RG079
      90,  212;  % RG080
     151,  183;  % RG081
      81,  400;  % RG082
      66,  191;  % RG083
      60,  264;  % RG087
      47,  345;  % RG088
      36,  323;  % RG089
      53,  328;  % RG090
      54,  379;  % RG091
      79,   35;  % RG092
     220,  262;  % RG100
      69,  190;  % RG384
      93,  224;  % RG398
    ];

% Derived fields (always computed for backward compatibility)
rainfall.use_uniform = ~strcmp(rainfall.type, 'real_event');
rainfall.intensity_ms = rainfall.uniform_intensity / 1000 / 3600;  % mm/hr -> m/s

%% Infiltration Parameters
% Walnut Gulch: Semi-arid soils, Green-Ampt infiltration
infiltration.model = 'green_ampt';  % Green-Ampt infiltration model
infiltration.enabled = true;        % Include infiltration state (F_cum) in model

% Green-Ampt parameters - spatially variable from 26 SSURGO soil map units (IMPLEMENTED)
% Representative values (spatially variable in grid):
infiltration.Ks_mean = 1.24e-5;    % Mean saturated hydraulic conductivity (m/s)
infiltration.Ks_min = 3.67e-6;     % Minimum Ks (Loam) (m/s)
infiltration.Ks_max = 5.83e-5;     % Maximum Ks (Sand) (m/s)
infiltration.psi_mean = 0.089;     % Mean suction head (m)
infiltration.psi_min = 0.049;      % Minimum psi (Sand) (m)
infiltration.psi_max = 0.316;      % Maximum psi (Clay) (m)
infiltration.theta_i = 0.15;       % Initial moisture content (-) uniform for dry season
infiltration.theta_s_mean = 0.445; % Mean saturated moisture content (-)
infiltration.theta_s_min = 0.437;  % Minimum porosity (Sand, Loamy Sand)
infiltration.theta_s_max = 0.463;  % Maximum porosity (Loam)
infiltration.delta_theta_mean = infiltration.theta_s_mean - infiltration.theta_i;

% Spatial variability statistics (from SSURGO mapping)
% Mapped from 26 SSURGO soil map units
% Result: 15.9x spatial variation in Ks (3.67e-6 to 5.83e-5 m/s)
% Coverage: 100% of grid cells (291/291)
% Data file: data/walnut_gulch/processed/ssurgo_parameters.csv

infiltration.use_spatial = true;    % Grid has spatially variable Ks, psi, theta_sat
infiltration.spatial_file = domain.grid_file;  % Parameters stored in grid file

% Calibration factor for Ks (multiplier)
% SSURGO lab-measured Ks overestimates field effective Ks due to:
%   - Surface sealing/crusting during intense rain
%   - Desert pavement and biological soil crusts
%   - Rock fragment armoring reducing infiltrable area
% Published calibrated values for WG models (KINEROS2) use Ks ~2-10 mm/hr
% vs SSURGO mean ~45 mm/hr → factor ~0.05-0.2
% Set to 1.0 to use raw SSURGO values (no calibration)
infiltration.Ks_factor = 0.05;  % Calibrated: KINEROS2 uses Ks ~2-10 mm/hr vs SSURGO ~45 mm/hr

% Antecedent soil moisture — cumulative infiltration prior to event
% Aug 23 follows Aug 16 storm by 7 days; partial recovery in semi-arid soil
% F_initial > F_min makes early runoff easier (soil partially wet)
infiltration.F_initial = 0.015;  % 15 mm antecedent infiltration (moderate)

%% Time Configuration
% Adapt duration to rainfall mode
time.start = 0.0;             % Start time (hours)
if strcmp(rainfall.type, 'real_event')
    % Real event: rainfall + recession
    if isfinite(rainfall.event_duration_min)
        time.end = rainfall.event_duration_min / 60;  % Match requested window
    else
        time.end = 25.0;       % 25 hours — captures full recession
    end
    time.dt_output = 1/60;     % 1-minute output (matches event data resolution)
elseif strcmp(rainfall.type, 'uniform_schedule')
    % Schedule mode: run until end of last pulse + 50% recession
    last_t = max(rainfall.schedule(:, 2));  % Last pulse end time [hr]
    time.end = last_t + 0.5 * last_t;      % Add 50% for recession
    time.dt_output = 0.05;     % 3-minute output
else
    % uniform_constant: rain duration + recession buffer
    time.end = rainfall.duration + 0.75;  % Storm + 45 min recession
    time.dt_output = 0.05;     % 3-minute output
end

% Convert to seconds for solver
time.start_s = time.start * 3600;
time.end_s = time.end * 3600;
time.dt_output_s = time.dt_output * 3600;

% Solver tolerances (different for ODE vs DAE)
% Local Inertial ODE is more stable, so can use tighter tolerances
if exist('hydraulic', 'var') && strcmp(hydraulic.routing_method, 'local_inertial')
    % Local Inertial ODE: can use standard tolerances
    time.RelTol = 1e-5;       % Standard relative tolerance for ODE
    time.AbsTol = 1e-6;       % Standard absolute tolerance for ODE
    time.MaxStep = 10;        % Smaller time steps for momentum accuracy
else
    % Manning CA DAE: need very relaxed tolerances for steep terrain
    time.RelTol = 1e-3;       % Very relaxed relative tolerance
    time.AbsTol = 1e-4;       % Very relaxed absolute tolerance  
    time.MaxStep = 10;        % Larger time steps allowed for algebraic fluxes
end

%% Numerical regularization for extremely steep terrain
% Override defaults for aggressive numerical stabilization
hydraulic.h_smooth = 5e-3;        % Very aggressive depth regularization (5mm vs default 2mm)
hydraulic.h_flow_threshold = 2e-3; % Threshold for flow computation (2mm)

%% DAE Configuration
% For Manning CA: 6-state DAE [h, Q_l, Q_r, Q_u, Q_d, F]
% Define state components as parameters
dae.n_base_states = 1;           % h (water depth)
dae.n_algebraic_states = 4;      % Q_l, Q_r, Q_u, Q_d (Manning CA fluxes)
dae.n_infiltration_states = 1;   % F_cum (Green-Ampt infiltration)

% Set formulation based on routing method
if strcmp(hydraulic.routing_method, 'manning_ca')
    if infiltration.enabled
        dae.formulation = '6state';  % Manning CA with infiltration
    else
        dae.formulation = '5state';  % Manning CA without infiltration
    end
    % Number of states: h + 4 fluxes + infiltration if enabled
    dae.n_states_per_cell = dae.n_base_states + dae.n_algebraic_states + double(infiltration.enabled) * dae.n_infiltration_states;
else
    error('Walnut Gulch currently only supports manning_ca routing method');
end

dae.solver = 'ode15i';               % MATLAB implicit DAE solver
dae.use_consistent_ic = true;        % Compute consistent initial conditions

% Solver method for Manning CA formulation
% Options:
%   'bdf'        - BDF fixed-point iteration 
%   'bdf_newton' - Newton
%   'ode15s'     - Semi-explicit DAE with regularization (robust for steep terrain)
dae.solver_method = 'bdf';

% Regularization for ode15s: 0 = pure DAE (may fail), >0 = regularized ODE
% Small epsilon converts algebraic constraints to fast dynamics: eps*dQ/dt = Q_manning - Q
% Recommended: 0.01 (10ms time constant) for steep terrain
dae.epsilon = 0.01;  % 10ms time constant for Q dynamics

% BDF settings (used for all bdf* methods)
dae.bdf_dt = 1;             % Time step [s] - increased from 0.25 (Newton converges in 2 iter)
dae.bdf_order = 1;           % 1 = Backward Euler (more stable), 2 = BDF2 (more accurate)
dae.bdf_newton_tol = 1e-5;   % Convergence tolerance (practical for steep terrain)
dae.bdf_max_iter = 100;       % Max fixed-point iterations (needs more than Newton)
dae.use_chordal = false;     % Not used for fixed-point BDF

% DAE system size (computed from domain)
dae.n_cells = domain.nx * domain.ny;  % Updated when DEM loaded
dae.n_total = dae.n_states_per_cell * dae.n_cells;

%% Observation Network (Real Gauges)
% Walnut Gulch has 19 instrumented flumes across the watershed.
% When using real event data, all 18 event flumes are configured automatically.
% For synthetic runs, observations remain empty (0 stations) or manually set.

obs.type = 'discharge';       % Observe discharge at stream gauges
obs.source = 'USDA_ARS';      % USDA Agricultural Research Service

% Observation uncertainty (based on rating curves and measurement error)
obs.std_depth = 0.02;         % Standard deviation for depth (m)
obs.std_discharge = 0.10;     % Standard deviation for discharge (m³/s) - 10%
obs.frequency = 1/60;         % Observation frequency (hours) - 1 minute

if strcmp(rainfall.type, 'real_event')
    % ---- Real event: 5 flumes from Aug 23, 2022 storm (50m grid) ----
    % Mapped from USDA-ARS instrumentation shapefiles (NAD83 UTM Zone 12N)
    % Source: data/walnut_gulch/raw/Files_SWRC/instrumentation/flumes/
    obs.n_stations = 5;
    obs.gauge_names = {'FL001', 'FL002', 'FL006', 'FL009', 'FL011'};
    obs.gauge_event_nums = [1, 2, 6, 9, 11];
    obs.locations = [...  % [row, col] — 1-indexed, row 1=north (50m grid)
         100,    1;  % FL001 — main outlet (Flume 1)
          87,  104;  % FL002
         108,  185;  % FL006
         121,  243;  % FL009
          70,  301;  % FL011
        ];
    obs.data_file = 'data/walnut_gulch/processed/event_20220823/flume_hydrograph_20220823.csv';
    obs.use_real_data = true;
    obs.use_synthetic = false;
else
    % ---- Synthetic / testing mode ----
    obs.gauge_names = {};         % Populate manually if needed
    obs.n_stations = 0;
    obs.locations = [];
    obs.data_file = '';
    obs.use_real_data = false;
    obs.use_synthetic = true;
end

%% State Estimation Configuration
estimation.method = 'closed_form';  %  closed-form algorithm
estimation.use_covariance = true;   % Use process covariance matrix

% Process noise parameters (tuned for real-world application)
estimation.sigma_h = 0.005;   % Depth process noise std (m)
estimation.sigma_q = 0.05;    % Discharge process noise std (m³/s)
estimation.sigma_v = 0.05;    % Velocity process noise std (m/s)
estimation.sigma_F = 0.001;   % Infiltration process noise std (m)

% Measurement noise (matches observation uncertainty)
estimation.R_depth = obs.std_depth^2;
estimation.R_discharge = obs.std_discharge^2;

% Spatial correlation (if used)
estimation.use_spatial_correlation = false;  % Simplified for computational efficiency
estimation.correlation_length = 500;         % Correlation length scale (m)

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
output.directory = 'results/walnut_gulch/AWR2026';
output.save_states = true;     % Save full state variables
output.save_diagnostics = true; % Save solver diagnostics
output.save_observations = true; % Save observations (real or synthetic)
output.save_estimates = true;   % Save state estimates

% Variables to save
output.variables = {'depth', 'discharge_x', 'discharge_y', ...
                    'velocity_x', 'velocity_y', 'infiltration'};

% Spatial outputs
output.save_spatial_maps = true;  % Save 2D spatial maps
output.map_times = [0.5, 1.0, 1.5, 2.0, 3.0];  % Times to save maps (hours)

% Save at gauge locations
output.save_gauge_timeseries = true;  % Detailed timeseries at gauges

%% Visualization Settings
viz.create_plots = true;
viz.plot_study_area = true;  % Study area grid/DEM visualization
viz.plot_pse_uncertainty_bounds = true;  % PSE uncertainty bounds (time series)
viz.plot_pse_spatial_uncertainty = true; % PSE spatial uncertainty maps
viz.plot_hydrograph = true;     % Hydrographs at gauges
viz.plot_spatial = true;        % Spatial maps at key times
viz.plot_comparison = true;     % DAE vs observations
viz.plot_estimation = true;     % State estimation results
viz.plot_uncertainty = true;    % Uncertainty bounds

% Figure export settings
viz.export_format = 'png';      % Options: 'png', 'pdf', 'eps'
viz.export_dpi = 300;           % Resolution for raster formats
viz.export_to_latex = true;     % Export using LaTeX utilities for paper

% Overlay options
viz.overlay_gauges = true;      % Show gauge locations on maps
viz.overlay_watershed = true;   % Show watershed boundary
viz.overlay_channels = true;    % Highlight main channels

%% Validation Settings
validation.run_hydropol2d = false;   % Compare with original HydroPol2D
validation.compare_observations = true;  % Compare with observed data
validation.compute_metrics = true;    % Compute error metrics

% Metrics to compute (against observations)
validation.metrics = {'RMSE', 'NSE', 'MAE', 'PeakError', 'VolumeError', ...
                      'TimeToPeak', 'PeakTime'};

% Performance evaluation
validation.evaluation_gauges = obs.gauge_names;  % Evaluate at all gauges
validation.evaluation_period = [time.start, time.end];  % Full simulation period

%% Data Loading Configuration
data.auto_load = true;          % Automatically load all required data
data.check_missing = true;      % Check for missing data and report

% Required data files (existing processed data)
data.required_files = {
    domain.spatial_data_file;   % Processed spatial data (DEM, bounds, etc.)
};

% Grid file (created if doesn't exist)
data.grid_file = domain.grid_file;
data.create_grid_if_missing = true;  % Run create_walnut_gulch_grid.m if needed

% Optional data files
data.optional_files = {
    'data/walnut_gulch/processed/WG_DAE_validation.mat';
    'data/walnut_gulch/processed/WG_DAE_PSE_results.mat';
};

%% Computational Settings
compute.use_parallel = false;   % Use parallel computing (if available)
compute.max_workers = 4;        % Maximum parallel workers

% Memory management
compute.save_incremental = true;  % Save results incrementally during simulation
compute.checkpoint_interval = 0.5;  % Checkpoint every 0.5 hours

%% Display Configuration Summary
fprintf('\n========================================\n');
fprintf('WALNUT GULCH DOMAIN CONFIGURATION\n');
fprintf('========================================\n');
fprintf('Watershed: %s (%s)\n', domain.name, domain.type);
fprintf('Location: Tombstone, Arizona, USA\n');
fprintf('Watershed area: %.1f km²\n', domain.area_km2);
fprintf('Grid resolution: %d x %d cells (%.0f m avg, CHECK GRID FILE!)\n', ...
    domain.nx, domain.ny, (domain.dx + domain.dy)/2);
fprintf('Spatial data: %s\n', domain.spatial_data_file);
fprintf('Grid file: %s\n', domain.grid_file);
fprintf('\n--- Spatially Variable Parameters ---\n');
if hydraulic.use_spatial_manning
    fprintf('Manning n: %.3f - %.3f (%.1fx variation, 24 veg types)\n', ...
        hydraulic.manning_n_min, hydraulic.manning_n_max, ...
        hydraulic.manning_n_max/hydraulic.manning_n_min);
else
    fprintf('Manning n: %.3f (uniform)\n', hydraulic.manning_n_default);
end
if infiltration.use_spatial
    fprintf('Hydraulic conductivity: %.2e - %.2e m/s (%.1fx variation, 26 soil units)\n', ...
        infiltration.Ks_min, infiltration.Ks_max, infiltration.Ks_max/infiltration.Ks_min);
    fprintf('Porosity: %.3f - %.3f\n', infiltration.theta_s_min, infiltration.theta_s_max);
    fprintf('Wetting front suction: %.3f - %.3f m\n', infiltration.psi_min, infiltration.psi_max);
else
    fprintf('Hydraulic conductivity: %.2e m/s (uniform)\n', infiltration.Ks_mean);
end
fprintf('\n--- Forcing and Simulation ---\n');
if strcmp(rainfall.type, 'real_event')
    fprintf('Rainfall: %s (%s)\n', rainfall.type, rainfall.event_description);
    fprintf('  Data file: %s\n', rainfall.event_data_file);
    fprintf('  Rain gages: %d mapped to grid (IDW: %s)\n', ...
        rainfall.n_gages, rainfall.interpolation);
    fprintf('  Mode: %s\n', ...
        ternary(rainfall.use_watershed_avg, 'watershed-average (Option A)', ...
        'spatial IDW (Option B)'));
else
    fprintf('Rainfall: %s (%.0f mm/hr, %.1f hr duration)\n', ...
        rainfall.type, rainfall.uniform_intensity, rainfall.duration_hours);
end
fprintf('Simulation time: %.1f to %.1f hours\n', time.start, time.end);
fprintf('DAE formulation: %s (%d states per cell)\n', ...
    dae.formulation, dae.n_states_per_cell);
fprintf('Total DAE system size: %d states\n', dae.n_total);
fprintf('\n--- Observation Network ---\n');
if obs.n_stations > 0
    fprintf('Observation stations: %d gauges\n', obs.n_stations);
    if obs.use_real_data
        fprintf('Using real observations: %s\n', obs.data_file);
    else
        fprintf('Using synthetic observations\n');
    end
else
    fprintf('Observation network: Not configured (19+ flumes available)\n');
    fprintf('  Set rainfall.type = ''real_event'' to auto-configure\n');
end
fprintf('========================================\n');

%% Helper function
function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
