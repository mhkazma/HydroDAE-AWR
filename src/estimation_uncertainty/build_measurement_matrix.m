function [C_local, K_vv_local, sensor_indices] = build_measurement_matrix(row, col, sensor_config)
% BUILD_MEASUREMENT_MATRIX  Construct local measurement matrix C for a cell
%
% USAGE:
%   [C_local, K_vv_local, sensor_indices] = build_measurement_matrix(row, col, sensor_config)
%
% DESCRIPTION:
%   For a given cell (row, col), constructs the local measurement matrix C
%   and measurement noise covariance K_vv following the PSE formulation
%   from :
%
%       y(k) = C * x(k) + v(k)
%
%   The measurement matrix C is a binary selection matrix that picks out
%   the measured states from the 6-state vector x = [h, Q_l, Q_r, Q_u, Q_d, F].
%
%   State-to-measurement mapping:
%     - Depth sensor  →  C row = [1, 0, 0, 0, 0, 0]  (selects h)
%     - Flow sensor (left)  →  C row = [0, 1, 0, 0, 0, 0]  (selects Q_l)
%     - Flow sensor (right) →  C row = [0, 0, 1, 0, 0, 0]  (selects Q_r)
%     - Flow sensor (up)    →  C row = [0, 0, 0, 1, 0, 0]  (selects Q_u)
%     - Flow sensor (down)  →  C row = [0, 0, 0, 0, 1, 0]  (selects Q_d)
%     - Infil sensor →  C row = [0, 0, 0, 0, 0, 1]  (selects F)
%
% INPUTS:
%   row           - Grid row index
%   col           - Grid column index
%   sensor_config - Sensor configuration from configure_sensors()
%
% OUTPUTS:
%   C_local       - n_y × 6 measurement matrix (n_y = sensors at this cell)
%                   Empty if no sensors at this cell
%   K_vv_local    - n_y × n_y measurement noise covariance (diagonal)
%   sensor_indices - n_y × 1 indices into state vector for each measurement
%
% REFERENCE:
%   Eq. (5): y(k) = C * x(k) + v(k)
%
% Author: M. H. Kazma
% Date: February 2026

n_states = 6;  % [h, Q_l, Q_r, Q_u, Q_d, F]

% Direction to state index mapping
% left=1→Q_l(2), right=2→Q_r(3), up=3→Q_u(4), down=4→Q_d(5)
dir_to_state = [2, 3, 4, 5];

% Count sensors at this cell
C_rows = [];
noise_vars = [];
s_indices = [];

% 1. Depth sensors at (row, col) → selects state 1 (h)
if sensor_config.n_depth > 0
    ds = sensor_config.depth_sensors;
    for si = 1:size(ds, 1)
        if ds(si, 1) == row && ds(si, 2) == col
            c_row = zeros(1, n_states);
            c_row(1) = 1;  % h
            C_rows = [C_rows; c_row]; %#ok<AGROW>
            noise_vars = [noise_vars; sensor_config.depth_noise_var(si)]; %#ok<AGROW>
            s_indices = [s_indices; 1]; %#ok<AGROW>
        end
    end
end

% 2. Flow sensors at (row, col, direction) → selects Q_l/Q_r/Q_u/Q_d
if sensor_config.n_flow > 0
    fs = sensor_config.flow_sensors;
    for si = 1:size(fs, 1)
        if fs(si, 1) == row && fs(si, 2) == col
            dir = fs(si, 3);
            state_idx = dir_to_state(dir);
            c_row = zeros(1, n_states);
            c_row(state_idx) = 1;
            C_rows = [C_rows; c_row]; %#ok<AGROW>
            noise_vars = [noise_vars; sensor_config.flow_noise_var(si)]; %#ok<AGROW>
            s_indices = [s_indices; state_idx]; %#ok<AGROW>
        end
    end
end

% 3. Infiltration sensors at (row, col) → selects state 6 (F)
if sensor_config.n_infil > 0
    is = sensor_config.infil_sensors;
    for si = 1:size(is, 1)
        if is(si, 1) == row && is(si, 2) == col
            c_row = zeros(1, n_states);
            c_row(6) = 1;  % F
            C_rows = [C_rows; c_row]; %#ok<AGROW>
            noise_vars = [noise_vars; sensor_config.infil_noise_var(si)]; %#ok<AGROW>
            s_indices = [s_indices; 6]; %#ok<AGROW>
        end
    end
end

% Assemble outputs
if isempty(C_rows)
    C_local = [];
    K_vv_local = [];
    sensor_indices = [];
else
    C_local = C_rows;
    K_vv_local = diag(noise_vars);
    sensor_indices = s_indices;
end

end
