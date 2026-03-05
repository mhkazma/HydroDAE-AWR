function [Z, params] = setup_vtilted_domain(dx, dy, nx_in, ny_in)
% setup_vtilted_domain: Create V-shaped tilted plane domain
%
% This function generates the topography and parameters for the V-Tilted
% benchmark case used in HydroPol2D validation studies.
%
% Inputs:
%   dx - Cell size in x-direction [m] (default: 20)
%   dy - Cell size in y-direction [m] (default: 20)
%   nx_in - Number of cells in x-direction (optional, overrides Lx/dx)
%   ny_in - Number of cells in y-direction (optional, overrides Ly/dy)
%
% Outputs:
%   Z - Elevation matrix [ny × nx] in meters
%   params - Structure containing domain parameters:
%       .grid_size - [ny, nx]
%       .dx, .dy - Cell dimensions [m]
%       .Z - Elevation [m]
%       .n_manning - Manning's n [s/m^(1/3)]
%       .g - Gravity [m/s²]
%       .h_min - Minimum water depth [m]
%       .h_0 - Initial water depth [m]
%       .outlet_indices - Indices of outlet cells
%       .outlet_type - 1 for free drainage
%       .outlet_direction - 'down', 'up', 'left', or 'right'
%       .slope_outlet - Slope at outlet boundary
%
% Example:
%   % Coarse resolution for testing
%   [Z, params] = setup_vtilted_domain(80, 80);
%
%   % Custom grid size
%   [Z, params] = setup_vtilted_domain(20, 20, 10, 10);
%   [Z, params] = setup_vtilted_domain(20, 20);
%
% Domain characteristics:
%   - Total domain: 1600m × 1000m
%   - V-shaped valley with channel in center
%   - Longitudinal slope: 2% (downward in y-direction)
%   - Lateral slope: 5% (toward center from both sides)
%   - Manning's n: 0.015 (hillslopes), 0.15 (channel)
%
% Author: M. H. Kazma
% Date: November 8, 2025

%% Parse inputs
if nargin < 1
    dx = 20;  % Default 20m resolution
end
if nargin < 2
    dy = 20;
end

%% Domain dimensions
if nargin >= 4 && ~isempty(nx_in) && ~isempty(ny_in)
    % Use provided grid dimensions
    nx = nx_in;
    ny = ny_in;
    Lx = nx * dx;
    Ly = ny * dy;
else
    % Use default domain size
    Lx = 1620;  % Domain length in x [m] (1620 so dx=20 gives odd nx=81)
    Ly = 1000;  % Domain length in y [m]
    nx = round(Lx / dx);
    ny = round(Ly / dy);
end

% Force odd lateral cell count for channel symmetry
% With odd nx, cell ceil(nx/2) is exactly the center -> channel sits on a cell
if mod(nx, 2) == 0
    nx = nx + 1;
    Lx = nx * dx;
    fprintf('  NOTE: nx adjusted to %d (odd) for channel symmetry\n', nx);
end

fprintf('Setting up V-Tilted domain:\n');
fprintf('  Resolution: dx=%.1fm, dy=%.1fm\n', dx, dy);
fprintf('  Grid: %d × %d = %d cells\n', ny, nx, ny*nx);
fprintf('  Domain size: %.1f × %.1f m\n\n', Lx, Ly);

%% Create V-shaped topography
Z = zeros(ny, nx);
center_x = Lx / 2;

for i = 1:ny
    for j = 1:nx
        % Cell center coordinates
        x = (j - 0.5) * dx;  % Use cell centers for symmetric geometry
        y = (i - 0.5) * dy;
        
        % Longitudinal slope (2% downward)
        z_y = 20 - 0.02 * y;
        
        % Lateral slope (5% toward center, V-shape)
        dist_from_center = abs(x - center_x);
        z_x = -0.05 * (Lx/2 - dist_from_center);  % Symmetric V-shape
        
        Z(i,j) = z_y + z_x;
    end
end

fprintf('Elevation range: %.1f to %.1f m\n\n', min(Z(:)), max(Z(:)));

%% Manning's roughness
n_manning = 0.015 * ones(ny, nx);  % Hillslopes

% Channel cells (3 cells wide in center)
% With odd nx, ceil(nx/2) is the exact center cell
channel_center = ceil(nx / 2);
channel_width = 1;  % Cells on each side of center
channel_cells = max(1, channel_center-channel_width) : min(nx, channel_center+channel_width);
n_manning(:, channel_cells) = 0.15;  % Higher roughness in channel

% Verify symmetry
Z_diff = abs(Z(:, 1:channel_center) - fliplr(Z(:, channel_center:end)));
max_asym = max(Z_diff(:));
fprintf('Channel symmetry check: max |Z_left - Z_right| = %.2e m\n', max_asym);

fprintf('Manning roughness:\n');
fprintf('  Hillslopes: n = 0.015\n');
fprintf('  Channel (center %d cells, cols %d-%d): n = 0.15\n\n', ...
    length(channel_cells), channel_cells(1), channel_cells(end));

%% Build parameters structure
params = struct();

% Grid parameters
params.grid_size = [ny, nx];
params.dx = dx;
params.dy = dy;

% Topography and roughness
params.Z = Z;
params.n_manning = n_manning;

% Physical constants
params.g = 9.81;  % Gravity [m/s²]
params.h_min = 1e-6;  % Minimum depth [m]
params.h_0 = 0;  % Initial depth [m]

% Outlet boundary (bottom row)
params.outlet_indices = (ny-1)*nx + (1:nx);
params.outlet_type = 1;  % Free drainage
params.outlet_direction = 'down';
params.slope_outlet = 0.02;  % Same as longitudinal slope

fprintf('Outlet configuration:\n');
fprintf('  Location: Bottom boundary (%d cells)\n', length(params.outlet_indices));
fprintf('  Type: Free drainage\n');
fprintf('  Slope: %.2f%%\n\n', params.slope_outlet*100);

end
