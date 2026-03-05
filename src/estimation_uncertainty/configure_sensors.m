function sensor_config = configure_sensors(sim_params, varargin)
% CONFIGURE_SENSORS  Define sensor placement for PSE with measurements
%
% USAGE:
%   sensor_config = configure_sensors(sim_params)
%   sensor_config = configure_sensors(sim_params, 'Name', Value, ...)
%
% DESCRIPTION:
%   Creates a sensor configuration structure that defines where sensors
%   (depth gauges, flow gauges, soil moisture probes) are placed on the
%   grid. This is used by run_pse_with_measurements to construct the
%   measurement matrix C from the PSE paper:
%
%       y(k) = C * x(k) + v(k)
%
%   where C is a binary selection matrix mapping states to measurements
%   and v(k) is measurement noise.
%
% STATE VECTOR (per cell): x = [h, Q_l, Q_r, Q_u, Q_d, F]
%   Index 1: h   - water depth [m]
%   Index 2: Q_l - left discharge [m^3/s]
%   Index 3: Q_r - right discharge [m^3/s]
%   Index 4: Q_u - up discharge [m^3/s]
%   Index 5: Q_d - down discharge [m^3/s]
%   Index 6: F   - cumulative infiltration [m]
%
% NAME-VALUE PAIRS:
%   'DepthSensors'     - Nx2 matrix of [row, col] for depth gauges
%   'DepthNoiseStd'    - Scalar or Nx1 noise std dev for depth [m]
%                        Default: 0.005 (5 mm)
%   'FlowSensors'      - Nx3 matrix of [row, col, direction]
%                        direction: 1=left, 2=right, 3=up, 4=down
%   'FlowNoiseStd'     - Scalar or Nx1 noise std dev for flow [m^3/s]
%                        Default: 1e-4
%   'InfilSensors'     - Nx2 matrix of [row, col] for infiltration probes
%   'InfilNoiseStd'    - Scalar or Nx1 noise std dev for F [m]
%                        Default: 0.002 (2 mm)
%   'OutletFlowSensor' - Logical, place flow sensor at outlet
%                        Default: true
%   'OutletNoiseStd'   - Scalar noise std dev for outlet flow [m^3/s]
%                        Default: 5e-4
%   'Strategy'         - Sensor placement strategy:
%                        'manual'  : Use provided sensor locations
%                        'uniform' : Place sensors on uniform grid
%                        'outlet'  : Only at outlet cell(s)
%                        'random'  : Random placement
%                        Default: 'manual'
%   'UniformSpacing'   - Spacing for 'uniform' strategy [cells]
%                        Default: 5
%   'NumRandom'        - Number of sensors for 'random' strategy
%                        Default: 10
%   'WeightMatrix'     - Diagonal weight matrix W for over-determined PSE
%                        Scalar or per-equation.
%                        Default: [] (unweighted)
%
% OUTPUT:
%   sensor_config - Structure with fields:
%     .depth_sensors    - Nx2 [row, col] locations
%     .depth_noise_var  - Nx1 noise variances (sigma^2)
%     .flow_sensors     - Mx3 [row, col, direction]
%     .flow_noise_var   - Mx1 noise variances
%     .infil_sensors    - Px2 [row, col] locations
%     .infil_noise_var  - Px1 noise variances
%     .n_depth          - Number of depth sensors
%     .n_flow           - Number of flow sensors
%     .n_infil          - Number of infiltration sensors
%     .n_total          - Total number of sensors
%     .weight_matrix    - Weight matrix W (or [])
%     .strategy         - Strategy used
%
% REFERENCE:
%   "Probabilistic State Estimation in Water Networks"
%   Equation (5): y(k) = C * x(k) + v(k)
%
% Author: M. H. Kazma
% Date: February 2026

%% Parse inputs
p = inputParser;
addRequired(p, 'sim_params', @isstruct);
addParameter(p, 'DepthSensors', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==2));
addParameter(p, 'DepthNoiseStd', 0.005, @isnumeric);
addParameter(p, 'FlowSensors', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==3));
addParameter(p, 'FlowNoiseStd', 1e-4, @isnumeric);
addParameter(p, 'InfilSensors', [], @(x) isempty(x) || (isnumeric(x) && size(x,2)==2));
addParameter(p, 'InfilNoiseStd', 0.002, @isnumeric);
addParameter(p, 'OutletFlowSensor', true, @islogical);
addParameter(p, 'OutletNoiseStd', 5e-4, @isnumeric);
addParameter(p, 'Strategy', 'manual', @(x) ismember(x, {'manual','uniform','outlet','random'}));
addParameter(p, 'UniformSpacing', 5, @isnumeric);
addParameter(p, 'NumRandom', 10, @isnumeric);
addParameter(p, 'WeightMatrix', [], @isnumeric);
parse(p, sim_params, varargin{:});
opts = p.Results;

ny = sim_params.grid_size(1);
nx = sim_params.grid_size(2);

% Active cell mask
if isfield(sim_params, 'use_compressed_indexing') && sim_params.use_compressed_indexing
    active_mask = false(ny, nx);
    ctg = sim_params.compressed_to_grid;
    active_mask(sub2ind([ny, nx], ctg(:,1), ctg(:,2))) = true;
else
    active_mask = true(ny, nx);
end

%% Generate sensor locations based on strategy
switch lower(opts.Strategy)
    case 'manual'
        depth_locs = opts.DepthSensors;
        flow_locs = opts.FlowSensors;
        infil_locs = opts.InfilSensors;
        
    case 'uniform'
        spacing = opts.UniformSpacing;
        rows_u = spacing:spacing:ny;
        cols_u = spacing:spacing:nx;
        [cc, rr] = meshgrid(cols_u, rows_u);
        depth_locs = [rr(:), cc(:)];
        % Filter to active cells
        valid = arrayfun(@(i) active_mask(depth_locs(i,1), depth_locs(i,2)), ...
            1:size(depth_locs,1));
        depth_locs = depth_locs(valid, :);
        flow_locs = opts.FlowSensors;
        infil_locs = opts.InfilSensors;
        
    case 'outlet'
        depth_locs = [];
        flow_locs = [];
        infil_locs = [];
        % Outlet sensor added below
        
    case 'random'
        n_rand = min(opts.NumRandom, sum(active_mask(:)));
        [ar, ac] = find(active_mask);
        perm = randperm(length(ar), n_rand);
        depth_locs = [ar(perm), ac(perm)];
        flow_locs = opts.FlowSensors;
        infil_locs = opts.InfilSensors;
        
    otherwise
        error('Unknown sensor strategy: %s', opts.Strategy);
end

%% Add outlet flow sensor if requested
if opts.OutletFlowSensor && isfield(sim_params, 'outlet_indices') && ~isempty(sim_params.outlet_indices)
    outlet_idx = sim_params.outlet_indices;
    if isfield(sim_params, 'outlet_direction')
        dir_str = sim_params.outlet_direction;
    else
        dir_str = 'right';
    end
    
    use_compressed = isfield(sim_params, 'use_compressed_indexing') && sim_params.use_compressed_indexing;
    
    % Direction mapping for fixed-direction modes
    dir_map = struct('left', 1, 'right', 2, 'up', 3, 'down', 4);
    
    if strcmp(dir_str, 'auto')
        % 'auto': detect outward-facing boundary faces PER CELL
        % A face is outward if the neighbor is outside the grid or inactive.
        % This creates one flow sensor per outward face per outlet cell,
        % matching the Q_outlet_ts computation in run_AWR2026_paper.m.
        n_faces_added = 0;
        
        for oi = 1:length(outlet_idx)
            if use_compressed
                ctg = sim_params.compressed_to_grid;
                gtc = sim_params.grid_to_compressed;
                o_row = ctg(outlet_idx(oi), 1);
                o_col = ctg(outlet_idx(oi), 2);
                
                % Left face (dir=1): neighbor col-1 is outside or inactive
                if o_col == 1 || gtc(o_row, o_col - 1) == 0
                    flow_locs = [flow_locs; o_row, o_col, 1]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
                % Right face (dir=2): neighbor col+1 is outside or inactive
                if o_col == nx || gtc(o_row, o_col + 1) == 0
                    flow_locs = [flow_locs; o_row, o_col, 2]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
                % Up face (dir=3): neighbor row-1 is outside or inactive
                if o_row == 1 || gtc(o_row - 1, o_col) == 0
                    flow_locs = [flow_locs; o_row, o_col, 3]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
                % Down face (dir=4): neighbor row+1 is outside or inactive
                if o_row == ny || gtc(o_row + 1, o_col) == 0
                    flow_locs = [flow_locs; o_row, o_col, 4]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
            else
                % Non-compressed grid: boundary = grid edge
                o_row = ceil(outlet_idx(oi) / nx);
                o_col = mod(outlet_idx(oi) - 1, nx) + 1;
                
                if o_col == 1
                    flow_locs = [flow_locs; o_row, o_col, 1]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
                if o_col == nx
                    flow_locs = [flow_locs; o_row, o_col, 2]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
                if o_row == 1
                    flow_locs = [flow_locs; o_row, o_col, 3]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
                if o_row == ny
                    flow_locs = [flow_locs; o_row, o_col, 4]; %#ok<AGROW>
                    n_faces_added = n_faces_added + 1;
                end
            end
        end
        fprintf('    auto outlet: %d cells, %d outward faces\n', length(outlet_idx), n_faces_added);
        
    elseif isfield(dir_map, dir_str)
        % Fixed direction mode: same direction for all outlet cells
        dir_num = dir_map.(dir_str);
        for oi = 1:length(outlet_idx)
            if use_compressed
                ctg = sim_params.compressed_to_grid;
                o_row = ctg(outlet_idx(oi), 1);
                o_col = ctg(outlet_idx(oi), 2);
            else
                o_row = ceil(outlet_idx(oi) / nx);
                o_col = mod(outlet_idx(oi) - 1, nx) + 1;
            end
            new_flow = [o_row, o_col, dir_num];
            if isempty(flow_locs)
                flow_locs = new_flow;
            else
                flow_locs = [flow_locs; new_flow]; %#ok<AGROW>
            end
        end
        fprintf('    fixed outlet direction -> %d (1=L,2=R,3=U,4=D)\n', dir_num);
    else
        warning('Unknown outlet_direction ''%s'', defaulting to down (4)', dir_str);
        dir_num = 4;
        for oi = 1:length(outlet_idx)
            if use_compressed
                ctg = sim_params.compressed_to_grid;
                o_row = ctg(outlet_idx(oi), 1);
                o_col = ctg(outlet_idx(oi), 2);
            else
                o_row = ceil(outlet_idx(oi) / nx);
                o_col = mod(outlet_idx(oi) - 1, nx) + 1;
            end
            flow_locs = [flow_locs; o_row, o_col, dir_num]; %#ok<AGROW>
        end
    end
end

%% Build noise variance vectors
if isempty(depth_locs)
    depth_noise_var = [];
    n_depth = 0;
else
    n_depth = size(depth_locs, 1);
    if isscalar(opts.DepthNoiseStd)
        depth_noise_var = (opts.DepthNoiseStd)^2 * ones(n_depth, 1);
    else
        depth_noise_var = opts.DepthNoiseStd(:).^2;
    end
end

if isempty(flow_locs)
    flow_noise_var = [];
    n_flow = 0;
else
    n_flow = size(flow_locs, 1);
    % Use outlet noise for outlet sensors, general noise for others
    if isscalar(opts.FlowNoiseStd)
        flow_noise_var = (opts.FlowNoiseStd)^2 * ones(n_flow, 1);
    else
        flow_noise_var = opts.FlowNoiseStd(:).^2;
    end
    % Override for outlet sensors
    if opts.OutletFlowSensor && n_flow > 0
        for fi = 1:n_flow
            is_at_outlet = false;
            if isfield(sim_params, 'outlet_indices')
                oi_row = flow_locs(fi, 1);
                oi_col = flow_locs(fi, 2);
                if use_compressed
                    cell_ci = sim_params.grid_to_compressed(oi_row, oi_col);
                else
                    cell_ci = (oi_row - 1) * nx + oi_col;
                end
                is_at_outlet = any(sim_params.outlet_indices == cell_ci);
            end
            if is_at_outlet
                flow_noise_var(fi) = (opts.OutletNoiseStd)^2;
            end
        end
    end
end

if isempty(infil_locs)
    infil_noise_var = [];
    n_infil = 0;
else
    n_infil = size(infil_locs, 1);
    if isscalar(opts.InfilNoiseStd)
        infil_noise_var = (opts.InfilNoiseStd)^2 * ones(n_infil, 1);
    else
        infil_noise_var = opts.InfilNoiseStd(:).^2;
    end
end

%% Assemble output structure
sensor_config = struct();
sensor_config.depth_sensors = depth_locs;
sensor_config.depth_noise_var = depth_noise_var;
sensor_config.flow_sensors = flow_locs;
sensor_config.flow_noise_var = flow_noise_var;
sensor_config.infil_sensors = infil_locs;
sensor_config.infil_noise_var = infil_noise_var;
sensor_config.n_depth = n_depth;
sensor_config.n_flow = n_flow;
sensor_config.n_infil = n_infil;
sensor_config.n_total = n_depth + n_flow + n_infil;
sensor_config.weight_matrix = opts.WeightMatrix;
sensor_config.strategy = opts.Strategy;
sensor_config.active_mask = active_mask;
sensor_config.grid_size = [ny, nx];

%% Validate sensor locations
validate_sensor_locations(sensor_config, active_mask, ny, nx);

fprintf('  Sensor configuration (%s strategy):\n', opts.Strategy);
fprintf('    Depth sensors:  %d\n', n_depth);
fprintf('    Flow sensors:   %d\n', n_flow);
fprintf('    Infil sensors:  %d\n', n_infil);
fprintf('    Total sensors:  %d\n', n_depth + n_flow + n_infil);
if n_depth > 0
    fprintf('    Depth noise std:  %.4f m (%.1f mm)\n', sqrt(depth_noise_var(1)), sqrt(depth_noise_var(1))*1000);
end
if n_flow > 0
    fprintf('    Flow noise std:   %.6f m^3/s\n', sqrt(flow_noise_var(1)));
end
if n_infil > 0
    fprintf('    Infil noise std:  %.4f m (%.1f mm)\n', sqrt(infil_noise_var(1)), sqrt(infil_noise_var(1))*1000);
end

end

function validate_sensor_locations(sc, active_mask, ny, nx)
% Check all sensor locations are within bounds and on active cells
locs_to_check = {};
names = {};
if sc.n_depth > 0
    locs_to_check{end+1} = sc.depth_sensors;
    names{end+1} = 'Depth';
end
if sc.n_flow > 0
    locs_to_check{end+1} = sc.flow_sensors(:, 1:2);
    names{end+1} = 'Flow';
end
if sc.n_infil > 0
    locs_to_check{end+1} = sc.infil_sensors;
    names{end+1} = 'Infil';
end

for li = 1:length(locs_to_check)
    locs = locs_to_check{li};
    for si = 1:size(locs, 1)
        r = locs(si, 1);
        c = locs(si, 2);
        if r < 1 || r > ny || c < 1 || c > nx
            error('%s sensor %d at (%d,%d) is out of grid bounds [%d x %d].', ...
                names{li}, si, r, c, ny, nx);
        end
        if ~active_mask(r, c)
            warning('%s sensor %d at (%d,%d) is on an inactive cell.', ...
                names{li}, si, r, c);
        end
    end
end
end
