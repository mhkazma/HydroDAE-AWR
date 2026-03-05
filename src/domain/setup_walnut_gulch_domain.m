function [Z, params] = setup_walnut_gulch_domain(grid_file)
% setup_walnut_gulch_domain: Load Walnut Gulch domain from grid file
%
% This function loads the pre-processed grid data for Walnut Gulch
% Experimental Watershed and prepares parameters for DAE simulation.
%
% Inputs:
%   grid_file - Path to grid .mat file (default: auto-detect from params)
%               Example: 'data/walnut_gulch/processed/walnut_gulch_grid_25x25.mat'
%
% Outputs:
%   Z - Elevation matrix [ny × nx] in meters
%   params - Structure containing domain parameters:
%       .grid_size - [ny, nx]
%       .dx, .dy - Cell dimensions [m]
%       .Z - Elevation [m]
%       .n_manning - Spatially variable Manning's n [s/m^(1/3)]
%       .Ks - Saturated hydraulic conductivity [m/s] (spatially variable)
%       .psi - Wetting front suction head [m] (spatially variable)
%       .theta_s - Saturated moisture content [-] (spatially variable)
%       .theta_i - Initial moisture content [-] (spatially variable)
%       .g - Gravity [m/s²]
%       .h_min - Minimum water depth [m]
%       .h_0 - Initial water depth [m]
%       .outlet_indices - Indices of outlet cells
%       .outlet_type - Boundary condition type
%       .valid_cells - Number of cells inside watershed
%       .X_grid, .Y_grid - Coordinate arrays
%       .slope_x, .slope_y, .slope_magnitude - Terrain slopes
%       .neighbors, .has_neighbor - Grid connectivity (4D)
%       .boundary_type - Boundary condition classification
%
% Example:
%   % Load default 25×25 grid
%   [Z, params] = setup_walnut_gulch_domain();
%
%   % Load specific grid file
%   grid_file = 'data/walnut_gulch/processed/walnut_gulch_grid_50x50.mat';
%   [Z, params] = setup_walnut_gulch_domain(grid_file);
%
% Notes:
%   - Grid file must be created by create_walnut_gulch_grid.m first
%   - Contains spatially variable parameters:
%     * Manning's n: 24 vegetation types (0.015 - 0.078)
%     * Ks: 26 SSURGO soil units (3.67e-6 - 5.83e-5 m/s)
%     * Green-Ampt parameters from literature
%   - Irregular domain with NaN values outside watershed boundary
%   - Outlet cells identified from lowest elevations at boundary
%
% Author: M. H. Kazma
% Date: November 11, 2025

%% Parse inputs
if nargin < 1 || isempty(grid_file)
    % Try to auto-detect grid file
    grid_file = 'data/walnut_gulch/processed/walnut_gulch_grid_25x25.mat';
end

%% Load grid data
fprintf('========================================\n');
fprintf('LOADING WALNUT GULCH DOMAIN\n');
fprintf('========================================\n');
fprintf('Grid file: %s\n', grid_file);

if ~exist(grid_file, 'file')
    error(['Grid file not found: %s\n' ...
           'Please run create_walnut_gulch_grid.m first to generate the grid.'], ...
           grid_file);
end

% Load all grid variables
grid_data = load(grid_file);

fprintf('Grid loaded successfully\n');
fprintf('  Grid size: %d × %d cells\n', grid_data.grid_size(1), grid_data.grid_size(2));
fprintf('  Valid cells: %d (%.1f%%)\n', grid_data.valid_cells, ...
    100*grid_data.valid_cells/prod(grid_data.grid_size));
fprintf('  Cell size: %.1f × %.1f m (avg: %.1f m)\n', ...
    grid_data.cell_size_x, grid_data.cell_size_y, grid_data.cell_size);

%% Extract elevation
Z = grid_data.Z_grid;

% Fill DEM pits/depressions to ensure drainage connectivity
% Pits cause unrealistic water accumulation (h >> 1m in interior cells)
valid_mask = ~isnan(Z);
fprintf('\nFilling DEM pits (priority flood)...\n');
Z = fill_dem_pits(Z, valid_mask);

% Get elevation statistics for valid cells
Z_valid = Z(~isnan(Z));
fprintf('\nElevation (after pit fill):\n');
fprintf('  Min: %.1f m\n', min(Z_valid));
fprintf('  Max: %.1f m\n', max(Z_valid));
fprintf('  Range: %.1f m\n', max(Z_valid) - min(Z_valid));

%% Extract slopes
slope_valid = grid_data.slope_magnitude(~isnan(grid_data.slope_magnitude));
fprintf('\nSlope:\n');
fprintf('  Mean: %.4f m/m (%.2f%%)\n', mean(slope_valid), mean(slope_valid)*100);
fprintf('  Max: %.4f m/m (%.2f%%)\n', max(slope_valid), max(slope_valid)*100);

%% Extract spatially variable parameters
% Manning's roughness (from 24 vegetation types)
n_manning = grid_data.n_manning;
n_valid = n_manning(~isnan(n_manning));
fprintf('\nManning Roughness (spatially variable):\n');
fprintf('  Range: %.3f - %.3f (%.1fx variation)\n', ...
    min(n_valid), max(n_valid), max(n_valid)/min(n_valid));
fprintf('  Mean: %.3f\n', mean(n_valid));
fprintf('  Source: 24 vegetation cover types\n');

% Hydraulic conductivity (from 26 SSURGO soil units)
Ks = grid_data.Ks;
Ks_valid = Ks(~isnan(Ks));
fprintf('\nHydraulic Conductivity (spatially variable):\n');
fprintf('  Range: %.2e - %.2e m/s (%.1fx variation)\n', ...
    min(Ks_valid), max(Ks_valid), max(Ks_valid)/min(Ks_valid));
fprintf('  Mean: %.2e m/s (%.1f mm/hr)\n', mean(Ks_valid), mean(Ks_valid)*3.6e6);
fprintf('  Source: 26 SSURGO soil map units\n');

% Wetting front suction
psi = grid_data.psi;
psi_valid = psi(~isnan(psi));
fprintf('\nWetting Front Suction:\n');
fprintf('  Range: %.3f - %.3f m\n', min(psi_valid), max(psi_valid));
fprintf('  Mean: %.3f m\n', mean(psi_valid));

% Saturated moisture content
theta_sat = grid_data.theta_sat;
theta_valid = theta_sat(~isnan(theta_sat));
fprintf('\nPorosity (theta_sat):\n');
fprintf('  Range: %.3f - %.3f\n', min(theta_valid), max(theta_valid));
fprintf('  Mean: %.3f\n', mean(theta_valid));

% Initial moisture content
theta_init = grid_data.theta_init;
theta_init_valid = theta_init(~isnan(theta_init));
if ~isempty(theta_init_valid)
    fprintf('  Initial: %.3f (uniform)\n', theta_init_valid(1));
else
    fprintf('  Initial: Not set\n');
end

% Effective porosity (eta = theta_s - theta_i)
eta = theta_sat - theta_init;
eta_valid = eta(~isnan(eta));
fprintf('\nEffective Porosity (eta = theta_s - theta_i):\n');
fprintf('  Range: %.3f - %.3f\n', min(eta_valid), max(eta_valid));
fprintf('  Mean: %.3f\n', mean(eta_valid));

fprintf('\nGreen-Ampt Parameters Summary:\n');
fprintf('  F_min: 10.0 mm (minimum infiltration for stability)\n');

%% Build parameters structure
params = struct();

% Grid parameters
params.grid_size = grid_data.grid_size;
ny = params.grid_size(1);
nx = params.grid_size(2);
params.dx = grid_data.cell_size_x;
params.dy = grid_data.cell_size_y;
params.cell_size = grid_data.cell_size;

% Topography
params.Z = Z;
params.slope_x = grid_data.slope_x;
params.slope_y = grid_data.slope_y;
params.slope_magnitude = grid_data.slope_magnitude;

% Coordinate arrays
params.X_grid = grid_data.X_grid;
params.Y_grid = grid_data.Y_grid;

% Spatially variable hydraulic parameters
params.n_manning = n_manning;

% Spatially variable infiltration parameters (Green-Ampt)
params.Ks = Ks;
params.psi = psi;
params.theta_s = theta_sat;
params.theta_i = theta_init;

% Derived infiltration parameters for DAE solver
% eta = effective porosity = theta_s - theta_i
params.eta = theta_sat - theta_init;  % Spatially variable
params.psi_f = psi;  % Alias for compatibility with add_infiltration_params
params.F_min = 0.010;  % Minimum cumulative infiltration [m] = 10mm (numerical floor)
params.F_initial = 0.010;  % Default: same as F_min (overridden by params file if set)

% Physical constants
params.g = 9.81;  % Gravity [m/s²]
params.h_min = 1e-6;  % Minimum depth [m]
params.h_0 = 0;  % Initial depth [m] - dry bed

% Slope statistics (for numerical regularization)
params.slope_max = max(slope_valid);
params.slope_mean = mean(slope_valid);

% Grid connectivity (4-directional: North, South, East, West)
params.neighbors = grid_data.neighbors;
params.has_neighbor = grid_data.has_neighbor;
params.connectivity = grid_data.connectivity;

% Boundary conditions
params.is_boundary = grid_data.is_boundary;
params.boundary_type = grid_data.boundary_type;
params.valid_cells = grid_data.valid_cells;
params.n_internal = grid_data.n_internal;
params.n_boundary = grid_data.n_boundary;
params.n_outlet = grid_data.n_outlet;

% Store outlet mask for later conversion (after grid_to_compressed is created)
outlet_mask = (grid_data.boundary_type == 2);  % Outlet cells (in full grid)

% Domain extent
params.watershed_bounds = grid_data.watershed_bounds;
params.watershed_area_km2 = grid_data.watershed_area_km2;

%% Create compressed indexing for valid cells only
% Build mapping: 2D grid (i,j) -> compressed valid cell index
% Only cells inside watershed boundary are included in state vector
valid_mask = ~isnan(Z);
params.valid_mask = valid_mask;

% Create mapping: grid_to_compressed(i,j) = compressed_index (0 if invalid)
params.grid_to_compressed = zeros(ny, nx);
compressed_idx = 0;
for i = 1:ny
    for j = 1:nx
        if valid_mask(i, j)
            compressed_idx = compressed_idx + 1;
            params.grid_to_compressed(i, j) = compressed_idx;
        end
    end
end

% Create inverse mapping: compressed_to_grid(k) = [i, j]
params.compressed_to_grid = zeros(params.valid_cells, 2);
compressed_idx = 0;
for i = 1:ny
    for j = 1:nx
        if valid_mask(i, j)
            compressed_idx = compressed_idx + 1;
            params.compressed_to_grid(compressed_idx, :) = [i, j];
        end
    end
end

fprintf('\nCompressed Indexing:\n');
fprintf('  Full grid cells: %d\n', ny*nx);
fprintf('  Valid cells only: %d\n', params.valid_cells);
fprintf('  Compression ratio: %.1f%%\n', 100*params.valid_cells/(ny*nx));

% Convert outlet indices to compressed space (now that grid_to_compressed exists)
outlet_indices_full = find(outlet_mask(:));
outlet_indices_compressed = [];
for idx = outlet_indices_full'
    [i, j] = ind2sub([ny, nx], idx);
    compressed_idx = params.grid_to_compressed(i, j);
    if compressed_idx > 0
        outlet_indices_compressed = [outlet_indices_compressed; compressed_idx];
    end
end

params.outlet_indices = outlet_indices_compressed;  % COMPRESSED indices
params.outlet_type = 1;  % Free drainage
params.outlet_direction = 'auto';  % Determined from flow accumulation

% Compute outlet slope from actual terrain slopes at outlet cells
if ~isempty(outlet_indices_compressed)
    outlet_slopes = zeros(length(outlet_indices_compressed), 1);
    for k = 1:length(outlet_indices_compressed)
        ci = outlet_indices_compressed(k);
        [oi, oj] = ind2sub([ny, nx], find(params.grid_to_compressed == ci));
        if ~isempty(oi)
            outlet_slopes(k) = params.slope_magnitude(oi(1), oj(1));
        end
    end
    % Use mean outlet slope, with minimum of 0.01 for numerical stability
    params.slope_outlet = max(0.01, mean(outlet_slopes(outlet_slopes > 0)));
else
    params.slope_outlet = 0.01;  % Default if no outlets
end
fprintf('  Outlet slope: %.4f\n', params.slope_outlet);

fprintf('  Outlet cells (full grid): %d\n', length(outlet_indices_full));
fprintf('  Outlet cells (compressed): %d\n', length(outlet_indices_compressed));

% DAE system dimensions (COMPRESSED - only valid cells)
% Manning CA: [h, Q_l, Q_r, Q_u, Q_d, F]
params.n_base_states = 1;           % h (water depth)
params.n_algebraic_states = 4;      % Q_l, Q_r, Q_u, Q_d (Manning CA fluxes)
params.n_infiltration_states = 1;   % F_cum (Green-Ampt infiltration)
params.n_states_per_cell = 6;       % Total: 1 + 4 + 1 = 6
params.n_states = params.valid_cells * params.n_states_per_cell;  % COMPRESSED SIZE
params.formulation = '6state';      % Manning CA with infiltration
params.routing_method = 'manning_ca';  % Routing method
params.use_compressed_indexing = true;  % Flag for residual function

%% Display boundary configuration
fprintf('\nBoundary Configuration:\n');
fprintf('  Internal cells: %d\n', params.n_internal);
fprintf('  Edge cells: %d\n', params.n_boundary - params.n_outlet);
fprintf('  Outlet cells: %d\n', params.n_outlet);
fprintf('  Outlet type: Free drainage\n');

%% Display DAE system size
fprintf('\nDAE System (Compressed):\n');
fprintf('  Formulation: %s\n', params.formulation);
fprintf('  Routing method: %s\n', params.routing_method);
fprintf('  States per cell: %d [h, Q_l, Q_r, Q_u, Q_d, F]\n', params.n_states_per_cell);
fprintf('    - Base states: %d (h)\n', params.n_base_states);
fprintf('    - Algebraic states: %d (Q_l, Q_r, Q_u, Q_d)\n', params.n_algebraic_states);
fprintf('    - Infiltration states: %d (F)\n', params.n_infiltration_states);
fprintf('  Active cells: %d (%.1f%% of grid)\n', params.valid_cells, 100*params.valid_cells/(ny*nx));
fprintf('  Total states: %d (was %d for full grid)\n', params.n_states, ny*nx*params.n_states_per_cell);
fprintf('  Grid formulation: 4-directional (Left, Right, Up, Down)\n');

fprintf('========================================\n');
fprintf('Domain setup complete\n');
fprintf('========================================\n\n');

end
