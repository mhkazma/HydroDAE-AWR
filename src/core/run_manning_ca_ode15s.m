function [t_out, x_out, results] = run_manning_ca_ode15s(params, R_input, duration, varargin)
% RUN_MANNING_CA_ODE15S Solve Manning CA as semi-explicit DAE using ode15s
%
% Uses the mass matrix formulation M*y' = f(t,y) where M is singular:
%   - Differential states (h, F): M diagonal = 1
%   - Algebraic states (Q_l, Q_r, Q_u, Q_d): M diagonal = 0
%
% This is a proper Index-1 semi-explicit DAE where the algebraic
% constraints enforce Q = Q_manning(h, slope) at each time step.
%
% Algebraic constraints:
%   0 = Q_manning(h, S) - Q
%
% State vector per cell: [h, Q_l, Q_r, Q_u, Q_d, F]
%   h    - Water depth [m] (differential)
%   Q_l  - Discharge to left neighbor [m^3/s] (algebraic)
%   Q_r  - Discharge to right neighbor [m^3/s] (algebraic)
%   Q_u  - Discharge to upper neighbor [m^3/s] (algebraic)
%   Q_d  - Discharge to lower neighbor [m^3/s] (algebraic)
%   F    - Cumulative infiltration [m] (differential)
%
% Inputs:
%   params    - Parameter structure from setup domain
%   R_input   - Rainfall intensity [m/s] (scalar or vector per cell)
%   duration  - Simulation duration [seconds]
%
% Optional name-value pairs:
%   'RelTol'          - Relative tolerance (default: 1e-4)
%   'AbsTol'          - Absolute tolerance (default: 1e-6)
%   'MaxStep'         - Maximum time step [s] (default: inf)
%   'DisplayProgress' - Show progress messages (default: true)
%
% Outputs:
%   t_out   - Time vector [s]
%   x_out   - State matrix (n_times x n_states)
%   results - Structure with water balance and diagnostics
%
% References:
%   Shampine, L.F., Reichelt, M.W., "The MATLAB ODE Suite", SIAM J. Sci.
%   Comput. 18 (1997), pp. 1-22.
%
% Author: M. H. Kazma
% Date: January 2026

%% Parse inputs
p = inputParser;
addParameter(p, 'RelTol', 1e-4);
addParameter(p, 'AbsTol', 1e-6);
addParameter(p, 'MaxStep', inf);
addParameter(p, 'DisplayProgress', true);
addParameter(p, 'Epsilon', 0);  % Regularization: 0 = pure DAE, >0 = singularly perturbed
parse(p, varargin{:});

RelTol = p.Results.RelTol;
AbsTol = p.Results.AbsTol;
MaxStep = p.Results.MaxStep;
display_progress = p.Results.DisplayProgress;
epsilon = p.Results.Epsilon;

%% Extract grid parameters
grid_size = params.grid_size;
nx = grid_size(2);
ny = grid_size(1);
dx = params.dx;
dy = params.dy;
cell_area = dx * dy;

% Check for compressed indexing (only valid cells)
if isfield(params, 'use_compressed_indexing') && params.use_compressed_indexing
    use_compressed = true;
    n_cells = params.valid_cells;
    grid_to_compressed = params.grid_to_compressed;
    compressed_to_grid = params.compressed_to_grid;
else
    use_compressed = false;
    n_cells = nx * ny;
end

% 6 states per cell for semi-explicit DAE
states_per_cell = 6;
n_states = n_cells * states_per_cell;

% Prepare rainfall vector
if isscalar(R_input)
    R = R_input * ones(n_cells, 1);
else
    R = R_input(:);
end

%% Display setup information
if display_progress
    fprintf('========================================\n');
    if epsilon == 0
        fprintf('Manning CA Semi-Explicit DAE (ode15s)\n');
        fprintf('========================================\n');
        fprintf('Formulation: M*y'' = f(t,y) with singular mass matrix\n');
    else
        fprintf('Manning CA Regularized ODE (ode15s)\n');
        fprintf('========================================\n');
        fprintf('Formulation: M*y'' = f(t,y) with eps=%.2e for Q\n', epsilon);
    end
    fprintf('Domain: %d x %d = %d cells\n', ny, nx, n_cells);
    fprintf('States: %d (6 per cell: h, Q_l, Q_r, Q_u, Q_d, F)\n', n_states);
    if epsilon == 0
        fprintf('  Differential: h, F (2 per cell)\n');
        fprintf('  Algebraic: Q_l, Q_r, Q_u, Q_d (4 per cell)\n');
    else
        fprintf('  Slow: h, F (2 per cell)\n');
        fprintf('  Fast: Q_l, Q_r, Q_u, Q_d (4 per cell, tau=%.2es)\n', epsilon);
    end
    fprintf('Duration: %.1f s (%.1f min)\n', duration, duration/60);
    fprintf('Tolerances: RelTol=%.0e, AbsTol=%.0e\n', RelTol, AbsTol);
    
    if isfield(params, 'Ks')
        if numel(params.Ks) > 1
            fprintf('Infiltration: Spatially variable Ks\n');
        else
            fprintf('Infiltration: Ks = %.2f mm/hr\n', params.Ks * 3600 * 1000);
        end
    end
    fprintf('\n');
end

%% Initialize state vector
x0 = zeros(n_states, 1);
F_min = params.F_min;
F_initial = F_min;
if isfield(params, 'F_initial')
    F_initial = max(params.F_initial, F_min);
end

for c = 1:n_cells
    base = (c-1) * states_per_cell;
    x0(base + 1) = 0;           % h = 0 (dry surface)
    x0(base + 2) = 0;           % Q_l = 0
    x0(base + 3) = 0;           % Q_r = 0
    x0(base + 4) = 0;           % Q_u = 0
    x0(base + 5) = 0;           % Q_d = 0
    x0(base + 6) = F_initial;   % F = F_initial (antecedent moisture)
end

if display_progress
    if F_initial > F_min
        fprintf('Initial: DRY surface, antecedent F=%.1f mm (F_min=%.1f mm)\n\n', F_initial*1000, F_min*1000);
    else
        fprintf('Initial conditions: h=0, Q=0, F=%.2f mm\n\n', F_min*1000);
    end
end

%% Build mass matrix
if display_progress
    fprintf('Building mass matrix...\n');
end

% M is diagonal: 1 for differential states (h, F)
% For algebraic states (Q): 0 for pure DAE, epsilon for regularization
M_diag = zeros(n_states, 1);
for c = 1:n_cells
    base = (c-1) * states_per_cell;
    M_diag(base + 1) = 1;        % h - differential
    M_diag(base + 2) = epsilon;  % Q_l - algebraic (or regularized)
    M_diag(base + 3) = epsilon;  % Q_r
    M_diag(base + 4) = epsilon;  % Q_u
    M_diag(base + 5) = epsilon;  % Q_d
    M_diag(base + 6) = 1;        % F - differential
end
M = spdiags(M_diag, 0, n_states, n_states);

n_diff = sum(M_diag == 1);
n_alg = sum(M_diag == 0);
n_fast = sum(M_diag > 0 & M_diag < 1);
if display_progress
    if epsilon == 0
        fprintf('  Differential equations: %d\n', n_diff);
        fprintf('  Algebraic equations: %d\n', n_alg);
    else
        fprintf('  Differential equations: %d\n', n_diff);
        fprintf('  Fast dynamics (eps=%.2e): %d\n', epsilon, n_fast);
    end
end

%% Build Jacobian sparsity pattern
if display_progress
    fprintf('Computing Jacobian sparsity pattern...\n');
end

Jpattern = compute_jacobian_pattern_semi_explicit(params);
nnz_J = nnz(Jpattern);
sparsity = 100 * (1 - nnz_J / n_states^2);

if display_progress
    fprintf('  Size: %d × %d\n', n_states, n_states);
    fprintf('  Non-zeros: %d (%.4f%% dense)\n', nnz_J, 100-sparsity);
end

%% Create RHS function handle
odefun = @(t, x) compute_rhs_semi_explicit(t, x, R, params);

%% Test RHS at initial conditions
if display_progress
    fprintf('\nTesting RHS at t=0...\n');
    f0 = odefun(0, x0);
    fprintf('  ||f(0,x0)|| = %.6e\n', norm(f0));
    fprintf('  max|f| = %.6e\n', max(abs(f0)));
    
    % Check consistency: algebraic residuals should be near zero
    alg_residual = 0;
    for c = 1:n_cells
        base = (c-1) * states_per_cell;
        alg_residual = max(alg_residual, abs(f0(base+2)));
        alg_residual = max(alg_residual, abs(f0(base+3)));
        alg_residual = max(alg_residual, abs(f0(base+4)));
        alg_residual = max(alg_residual, abs(f0(base+5)));
    end
    fprintf('  max|algebraic residual| = %.6e\n\n', alg_residual);
end

%% Configure ode15s options
% For DAEs, provide initial slope (y'(0)) to help the solver
% For consistent ICs: h'(0) = R (rainfall), Q(0) = 0, F'(0) = 0
yp0 = zeros(n_states, 1);
for c = 1:n_cells
    base = (c-1) * states_per_cell;
    yp0(base + 1) = R(c);  % dh/dt = R (initial, before infiltration starts)
    % Q derivatives are algebraic, set to 0
    yp0(base + 6) = 0;     % dF/dt = 0 (no infiltration at h=0)
end

% Set MassSingular based on epsilon
if epsilon == 0
    mass_singular = 'yes';
else
    mass_singular = 'no';  % Regularized system is not singular
end

opts = odeset(...
    'Mass', M, ...
    'MStateDependence', 'none', ...
    'MassSingular', mass_singular, ...
    'InitialSlope', yp0, ...
    'RelTol', RelTol, ...
    'AbsTol', AbsTol, ...
    'JPattern', Jpattern, ...
    'Stats', 'on');

if ~isinf(MaxStep)
    opts = odeset(opts, 'MaxStep', MaxStep);
    if display_progress
        fprintf('MaxStep: %.2f s\n', MaxStep);
    end
end

%% Solve the DAE
if display_progress
    fprintf('Calling ode15s...\n');
    fprintf('----------------------------------------\n');
end

tic;
try
    [t_out, x_out] = ode15s(odefun, [0, duration], x0, opts);
    comp_time = toc;
catch ME
    comp_time = toc;
    if display_progress
        fprintf('\n*** Solver failed after %.1f s ***\n', comp_time);
        fprintf('Error: %s\n', ME.message);
    end
    rethrow(ME);
end

%% Post-process results
if display_progress
    fprintf('----------------------------------------\n');
    fprintf('Simulation complete\n');
    fprintf('  Wall time: %.2f s\n', comp_time);
    fprintf('  Time steps: %d\n', length(t_out));
    fprintf('  Final time: %.1f s (%.1f min)\n\n', t_out(end), t_out(end)/60);
end

%% Extract final state and compute water balance
x_final = x_out(end, :)';

% Reconstruct spatial grids
h_final = zeros(ny, nx);
F_final = zeros(ny, nx);
Q_out_final = zeros(ny, nx);

if use_compressed
    for c = 1:n_cells
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
        base = (c-1) * states_per_cell;
        
        h_final(i, j) = max(x_final(base + 1), 0);
        F_final(i, j) = x_final(base + 6);
        
        % Total outflow from cell
        Q_out_final(i, j) = x_final(base + 2) + x_final(base + 3) + ...
                           x_final(base + 4) + x_final(base + 5);
    end
else
    for i = 1:ny
        for j = 1:nx
            c = (i-1)*nx + j;
            base = (c-1) * states_per_cell;
            
            h_final(i, j) = max(x_final(base + 1), 0);
            F_final(i, j) = x_final(base + 6);
            Q_out_final(i, j) = x_final(base + 2) + x_final(base + 3) + ...
                               x_final(base + 4) + x_final(base + 5);
        end
    end
end

%% Compute actual outlet flow through boundary
% Integrate discharge at outlet boundary over time
total_outflow = 0;
outlet_direction = 'none';
if isfield(params, 'outlet_direction')
    outlet_direction = params.outlet_direction;
end

% Get outlet cell discharges based on direction
for k = 1:length(t_out)
    x_k = x_out(k, :)';
    Q_outlet_k = 0;
    
    if use_compressed
        for c = 1:n_cells
            i = compressed_to_grid(c, 1);
            j = compressed_to_grid(c, 2);
            base = (c-1) * states_per_cell;
            
            % Check if this is an outlet boundary cell
            is_outlet = false;
            Q_boundary = 0;
            
            if strcmp(outlet_direction, 'down') && i == ny
                is_outlet = true;
                Q_boundary = x_k(base + 5);  % Q_d (downward)
            elseif strcmp(outlet_direction, 'up') && i == 1
                is_outlet = true;
                Q_boundary = x_k(base + 4);  % Q_u (upward)
            elseif strcmp(outlet_direction, 'left') && j == 1
                is_outlet = true;
                Q_boundary = x_k(base + 2);  % Q_l (leftward)
            elseif strcmp(outlet_direction, 'right') && j == nx
                is_outlet = true;
                Q_boundary = x_k(base + 3);  % Q_r (rightward)
            end
            
            if is_outlet
                Q_outlet_k = Q_outlet_k + max(Q_boundary, 0);
            end
        end
    else
        for i = 1:ny
            for j = 1:nx
                c = (i-1)*nx + j;
                base = (c-1) * states_per_cell;
                
                is_outlet = false;
                Q_boundary = 0;
                
                if strcmp(outlet_direction, 'down') && i == ny
                    is_outlet = true;
                    Q_boundary = x_k(base + 5);
                elseif strcmp(outlet_direction, 'up') && i == 1
                    is_outlet = true;
                    Q_boundary = x_k(base + 4);
                elseif strcmp(outlet_direction, 'left') && j == 1
                    is_outlet = true;
                    Q_boundary = x_k(base + 2);
                elseif strcmp(outlet_direction, 'right') && j == nx
                    is_outlet = true;
                    Q_boundary = x_k(base + 3);
                end
                
                if is_outlet
                    Q_outlet_k = Q_outlet_k + max(Q_boundary, 0);
                end
            end
        end
    end
    
    % Integrate using trapezoidal rule
    if k > 1
        dt = t_out(k) - t_out(k-1);
        total_outflow = total_outflow + 0.5 * (Q_outlet_k + Q_outlet_prev) * dt;
    end
    Q_outlet_prev = Q_outlet_k;
end

%% Compute water balance
total_rainfall = mean(R) * duration * n_cells * cell_area;

% Sum only over valid cells (where simulation happened)
if use_compressed
    valid_mask = params.valid_mask;
    h_valid = h_final(valid_mask);
    F_valid = F_final(valid_mask);
    total_storage = sum(h_valid) * cell_area;
    total_infiltration = sum(F_valid - F_initial) * cell_area;
else
    total_storage = sum(h_final(:)) * cell_area;
    total_infiltration = sum(F_final(:) - F_initial) * cell_area;
end

balance_error = abs(total_rainfall - total_storage - total_infiltration - total_outflow);
balance_error_pct = 100 * balance_error / max(total_rainfall, 1e-10);

if display_progress
    fprintf('Water Balance:\n');
    fprintf('  Rainfall:      %12.4f m^3 (100.0%%)\n', total_rainfall);
    fprintf('  Ponding:       %12.4f m^3 (%5.1f%%)\n', total_storage, 100*total_storage/max(total_rainfall,1e-10));
    fprintf('  Infiltration:  %12.4f m^3 (%5.1f%%)\n', total_infiltration, 100*total_infiltration/max(total_rainfall,1e-10));
    fprintf('  Runoff:        %12.4f m^3 (%5.1f%%)\n', total_outflow, 100*total_outflow/max(total_rainfall,1e-10));
    fprintf('  --------------------------------\n');
    fprintf('  Balance error: %12.4f m^3 (%5.2f%%)\n\n', balance_error, balance_error_pct);
end

%% Package results
results = struct();
results.t = t_out;
results.x = x_out;
results.h_final = h_final;
results.F_final = F_final;
results.Q_out_final = Q_out_final;
results.computation_time = comp_time;
results.comp_time = comp_time;  % Alias for backward compatibility
results.n_steps = length(t_out);

results.water_balance.rainfall = total_rainfall;
results.water_balance.ponding = total_storage;
results.water_balance.infiltration = total_infiltration;
results.water_balance.runoff = total_outflow;
results.water_balance.error_pct = balance_error_pct;

results.solver_info.method = 'ode15s';
results.solver_info.formulation = 'semi-explicit DAE';
results.solver_info.n_states = n_states;
results.solver_info.n_differential = n_diff;
results.solver_info.n_algebraic = n_alg;
results.solver_info.nnz_jacobian = nnz_J;

end


%% ========================================================================
%  RHS FUNCTION: Semi-explicit DAE formulation
%  ========================================================================
function f = compute_rhs_semi_explicit(t, x, R, params)
% COMPUTE_RHS_SEMI_EXPLICIT Compute f(t,y) for M*y' = f(t,y)
%
% For differential states (h, F):
%   f_h = dh/dt = R - f_infil + (Q_in - Q_out) / A
%   f_F = dF/dt = f_infil
%
% For algebraic states (Q):
%   f_Q = 0 = Q - Q_manning(h, S, ...)
%
% This is the semi-explicit form where differential equations are explicit
% and algebraic constraints are implicit.

%% Extract parameters
grid_size = params.grid_size;
nx = grid_size(2);
ny = grid_size(1);
dx = params.dx;
dy = params.dy;
cell_area = dx * dy;

Z = params.Z;
n_manning = params.n_manning;
h_0 = params.h_0;
F_min = params.F_min;

% Regularization parameters
h_smooth = 2e-3;      % Depth smoothing [m]
h_flow_min = 1e-4;    % Minimum flow depth [m]
v_max = 10.0;         % Maximum velocity [m/s]

% Infiltration parameters
if isfield(params, 'Ks')
    Ks = params.Ks;
    eta = params.eta;
    psi_f = params.psi_f;
    use_infiltration = true;
else
    use_infiltration = false;
    Ks = 0; eta = 0; psi_f = 0;
end

% Compressed indexing
if isfield(params, 'use_compressed_indexing') && params.use_compressed_indexing
    use_compressed = true;
    n_cells = params.valid_cells;
    grid_to_compressed = params.grid_to_compressed;
    compressed_to_grid = params.compressed_to_grid;
else
    use_compressed = false;
    n_cells = nx * ny;
    grid_to_compressed = [];
    compressed_to_grid = [];
end

% Outlet boundary parameters
if isfield(params, 'outlet_direction')
    outlet_direction = params.outlet_direction;
    slope_outlet = params.slope_outlet;
else
    outlet_direction = 'none';
    slope_outlet = 0;
end

states_per_cell = 6;
n_states = n_cells * states_per_cell;
f = zeros(n_states, 1);

%% Compute RHS for each cell
if use_compressed
    for c = 1:n_cells
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
        base = (c-1) * states_per_cell;
        
        % Extract states
        h = max(x(base + 1), 0);
        Q_l = x(base + 2);
        Q_r = x(base + 3);
        Q_u = x(base + 4);
        Q_d = x(base + 5);
        F = max(x(base + 6), F_min);
        
        % Local parameters
        z = Z(i, j);
        n = n_manning(i, j);
        
        % Effective depth with regularization
        h_reg = sqrt(h^2 + h_smooth^2);
        h_eff = max(h_reg - h_0, h_smooth);
        
        % Water surface elevation
        wse = z + h;
        
        % Get neighbor water surface elevations (with outlet BC)
        [wse_l, wse_r, wse_u, wse_d] = get_neighbor_wse_with_outlet(i, j, x, Z, ...
            grid_to_compressed, ny, nx, states_per_cell, wse, dx, dy, ...
            outlet_direction, slope_outlet);
        
        % Compute slopes (positive = downhill from this cell)
        S_l = (wse - wse_l) / dx;
        S_r = (wse - wse_r) / dx;
        S_u = (wse - wse_u) / dy;
        S_d = (wse - wse_d) / dy;
        
        % Compute Manning discharge for each direction
        Q_manning_l = compute_manning_Q(h_eff, S_l, n, dy, h, h_flow_min, v_max);
        Q_manning_r = compute_manning_Q(h_eff, S_r, n, dy, h, h_flow_min, v_max);
        Q_manning_u = compute_manning_Q(h_eff, S_u, n, dx, h, h_flow_min, v_max);
        Q_manning_d = compute_manning_Q(h_eff, S_d, n, dx, h, h_flow_min, v_max);
        
        % Compute inflow from neighbors
        Q_in = compute_inflow(i, j, x, grid_to_compressed, ny, nx, states_per_cell);
        
        % Total outflow
        Q_out = Q_l + Q_r + Q_u + Q_d;
        
        % Infiltration (Green-Ampt)
        if use_infiltration && h > h_flow_min
            if isscalar(Ks)
                Ks_local = Ks;
                psi_local = psi_f;
                eta_local = eta;
            else
                Ks_local = Ks(i, j);
                psi_local = psi_f(i, j);
                eta_local = eta(i, j);
            end
            
            % Wetting front potential
            wf = 1 + (psi_local * eta_local) / max(F, F_min);
            f_infil = Ks_local * wf;
            
            % Limit infiltration to available water
            f_infil = min(f_infil, h / 0.01);  % dt approximation
        else
            f_infil = 0;
        end
        
        % === Differential equations ===
        % dh/dt = R - f_infil + (Q_in - Q_out) / A
        f(base + 1) = R(c) - f_infil + (Q_in - Q_out) / cell_area;
        
        % === Algebraic equations (constraints) ===
        % For semi-explicit DAE: 0 = g(y) where g = Q_manning - Q
        % ode15s solves 0*Q' = Q_manning - Q, i.e., Q = Q_manning
        f(base + 2) = Q_manning_l - Q_l;
        f(base + 3) = Q_manning_r - Q_r;
        f(base + 4) = Q_manning_u - Q_u;
        f(base + 5) = Q_manning_d - Q_d;
        
        % === Differential equation ===
        % dF/dt = f_infil
        f(base + 6) = f_infil;
    end
else
    % Full grid (no compression)
    for i = 1:ny
        for j = 1:nx
            c = (i-1)*nx + j;
            base = (c-1) * states_per_cell;
            
            % Extract states
            h = max(x(base + 1), 0);
            Q_l = x(base + 2);
            Q_r = x(base + 3);
            Q_u = x(base + 4);
            Q_d = x(base + 5);
            F = max(x(base + 6), F_min);
            
            z = Z(i, j);
            n = n_manning(i, j);
            
            h_reg = sqrt(h^2 + h_smooth^2);
            h_eff = max(h_reg - h_0, h_smooth);
            wse = z + h;
            
            % Neighbor WSE with outlet boundary conditions
            wse_l = wse; wse_r = wse; wse_u = wse; wse_d = wse;
            
            if j > 1
                nc = (i-1)*nx + (j-1);
                h_n = max(x((nc-1)*states_per_cell + 1), 0);
                wse_l = Z(i, j-1) + h_n;
            elseif strcmp(outlet_direction, 'left')
                % Free drainage to left: fictitious neighbor at lower elevation
                wse_l = wse - slope_outlet * dx;
            end
            if j < nx
                nc = (i-1)*nx + (j+1);
                h_n = max(x((nc-1)*states_per_cell + 1), 0);
                wse_r = Z(i, j+1) + h_n;
            elseif strcmp(outlet_direction, 'right')
                wse_r = wse - slope_outlet * dx;
            end
            if i > 1
                nc = (i-2)*nx + j;
                h_n = max(x((nc-1)*states_per_cell + 1), 0);
                wse_u = Z(i-1, j) + h_n;
            elseif strcmp(outlet_direction, 'up')
                wse_u = wse - slope_outlet * dy;
            end
            if i < ny
                nc = i*nx + j;
                h_n = max(x((nc-1)*states_per_cell + 1), 0);
                wse_d = Z(i+1, j) + h_n;
            elseif strcmp(outlet_direction, 'down')
                % Free drainage downward: fictitious neighbor at lower elevation
                wse_d = wse - slope_outlet * dy;
            end
            
            % Slopes
            S_l = (wse - wse_l) / dx;
            S_r = (wse - wse_r) / dx;
            S_u = (wse - wse_u) / dy;
            S_d = (wse - wse_d) / dy;
            
            % Manning discharge
            Q_manning_l = compute_manning_Q(h_eff, S_l, n, dy, h, h_flow_min, v_max);
            Q_manning_r = compute_manning_Q(h_eff, S_r, n, dy, h, h_flow_min, v_max);
            Q_manning_u = compute_manning_Q(h_eff, S_u, n, dx, h, h_flow_min, v_max);
            Q_manning_d = compute_manning_Q(h_eff, S_d, n, dx, h, h_flow_min, v_max);
            
            % Inflow from neighbors
            Q_in = 0;
            if j > 1  % Left neighbor's right outflow
                nc = (i-1)*nx + (j-1);
                Q_in = Q_in + max(x((nc-1)*states_per_cell + 3), 0);
            end
            if j < nx  % Right neighbor's left outflow
                nc = (i-1)*nx + (j+1);
                Q_in = Q_in + max(x((nc-1)*states_per_cell + 2), 0);
            end
            if i > 1  % Upper neighbor's down outflow
                nc = (i-2)*nx + j;
                Q_in = Q_in + max(x((nc-1)*states_per_cell + 5), 0);
            end
            if i < ny  % Lower neighbor's up outflow
                nc = i*nx + j;
                Q_in = Q_in + max(x((nc-1)*states_per_cell + 4), 0);
            end
            
            Q_out = Q_l + Q_r + Q_u + Q_d;
            
            % Infiltration
            if use_infiltration && h > h_flow_min
                if isscalar(Ks)
                    Ks_local = Ks; psi_local = psi_f; eta_local = eta;
                else
                    Ks_local = Ks(i,j); psi_local = psi_f(i,j); eta_local = eta(i,j);
                end
                wf = 1 + (psi_local * eta_local) / max(F, F_min);
                f_infil = min(Ks_local * wf, h / 0.01);
            else
                f_infil = 0;
            end
            
            % Differential: dh/dt
            f(base + 1) = R(c) - f_infil + (Q_in - Q_out) / cell_area;
            
            % Algebraic: 0 = Q_manning - Q
            f(base + 2) = Q_manning_l - Q_l;
            f(base + 3) = Q_manning_r - Q_r;
            f(base + 4) = Q_manning_u - Q_u;
            f(base + 5) = Q_manning_d - Q_d;
            
            % Differential: dF/dt
            f(base + 6) = f_infil;
        end
    end
end

end


%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================

function Q = compute_manning_Q(h_eff, S, n, width, h, h_min, v_max)
% COMPUTE_MANNING_Q Compute Manning discharge with smooth regularization
%
% Uses smooth transitions to avoid discontinuities that cause solver failures.
% Critical for steep terrain where adaptive solvers need smooth Jacobians.

    % Smooth slope regularization: sqrt(S^2 + eps) instead of hard cutoff at S=0
    S_smooth = 1e-8;  % Small regularization to avoid zero slope singularity
    S_eff = sqrt(S^2 + S_smooth^2);  % Always positive, smooth at S=0
    
    % Only compute for positive (downhill) slope
    if S > 0
        % Velocity from Manning equation with regularized slope
        v = (1/n) * h_eff^(2/3) * sqrt(S_eff);
        v = min(v, v_max);
        
        % Smooth flow factor (0 to 1) - cubic smoothstep
        % Transitions smoothly from 0 at h=0 to 1 at h=10*h_min
        ff = min(1.0, max(0.0, h / (10 * h_min)));
        ff = ff * ff * (3 - 2*ff);  % Smoothstep: C1 continuous
        
        Q = ff * v * h_eff * width;
    else
        % No flow uphill - but use smooth transition near S=0
        Q = 0;
    end
end


function [wse_l, wse_r, wse_u, wse_d] = get_neighbor_wse_with_outlet(i, j, x, Z, ...
    grid_to_compressed, ny, nx, spc, wse_self, dx, dy, outlet_direction, slope_outlet)
% GET_NEIGHBOR_WSE_WITH_OUTLET Get water surface elevations with outlet BC
%
% For interior neighbors: uses actual water surface elevation
% For boundary cells at outlet: uses free drainage BC (fictitious lower WSE)
% For non-outlet boundaries: uses no-flow BC (same WSE)

    wse_l = wse_self;
    wse_r = wse_self;
    wse_u = wse_self;
    wse_d = wse_self;
    
    % Left neighbor
    if j > 1
        nc = grid_to_compressed(i, j-1);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_l = Z(i, j-1) + h_n;
        end
    elseif strcmp(outlet_direction, 'left')
        % Free drainage to left
        wse_l = wse_self - slope_outlet * dx;
    end
    
    % Right neighbor
    if j < nx
        nc = grid_to_compressed(i, j+1);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_r = Z(i, j+1) + h_n;
        end
    elseif strcmp(outlet_direction, 'right')
        wse_r = wse_self - slope_outlet * dx;
    end
    
    % Upper neighbor
    if i > 1
        nc = grid_to_compressed(i-1, j);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_u = Z(i-1, j) + h_n;
        end
    elseif strcmp(outlet_direction, 'up')
        wse_u = wse_self - slope_outlet * dy;
    end
    
    % Lower neighbor (i = ny is bottom boundary)
    if i < ny
        nc = grid_to_compressed(i+1, j);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_d = Z(i+1, j) + h_n;
        end
    elseif strcmp(outlet_direction, 'down')
        % Free drainage downward: fictitious WSE at outlet slope below
        wse_d = wse_self - slope_outlet * dy;
    end
end


function [wse_l, wse_r, wse_u, wse_d] = get_neighbor_wse(i, j, x, Z, ...
    grid_to_compressed, ny, nx, spc, wse_self)
% GET_NEIGHBOR_WSE Get water surface elevations of neighbors
%
% Returns wse_self for boundary cells (no-flow BC)

    wse_l = wse_self;
    wse_r = wse_self;
    wse_u = wse_self;
    wse_d = wse_self;
    
    % Left neighbor
    if j > 1
        nc = grid_to_compressed(i, j-1);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_l = Z(i, j-1) + h_n;
        end
    end
    
    % Right neighbor
    if j < nx
        nc = grid_to_compressed(i, j+1);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_r = Z(i, j+1) + h_n;
        end
    end
    
    % Upper neighbor
    if i > 1
        nc = grid_to_compressed(i-1, j);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_u = Z(i-1, j) + h_n;
        end
    end
    
    % Lower neighbor
    if i < ny
        nc = grid_to_compressed(i+1, j);
        if nc > 0
            h_n = max(x((nc-1)*spc + 1), 0);
            wse_d = Z(i+1, j) + h_n;
        end
    end
end


function Q_in = compute_inflow(i, j, x, grid_to_compressed, ny, nx, spc)
% COMPUTE_INFLOW Sum inflow from neighboring cells
%
% Inflow comes from:
%   - Left neighbor's Q_r (index 3)
%   - Right neighbor's Q_l (index 2)
%   - Upper neighbor's Q_d (index 5)
%   - Lower neighbor's Q_u (index 4)

    Q_in = 0;
    
    % From left neighbor (their Q_r)
    if j > 1
        nc = grid_to_compressed(i, j-1);
        if nc > 0
            Q_in = Q_in + max(x((nc-1)*spc + 3), 0);
        end
    end
    
    % From right neighbor (their Q_l)
    if j < nx
        nc = grid_to_compressed(i, j+1);
        if nc > 0
            Q_in = Q_in + max(x((nc-1)*spc + 2), 0);
        end
    end
    
    % From upper neighbor (their Q_d)
    if i > 1
        nc = grid_to_compressed(i-1, j);
        if nc > 0
            Q_in = Q_in + max(x((nc-1)*spc + 5), 0);
        end
    end
    
    % From lower neighbor (their Q_u)
    if i < ny
        nc = grid_to_compressed(i+1, j);
        if nc > 0
            Q_in = Q_in + max(x((nc-1)*spc + 4), 0);
        end
    end
end


%% ========================================================================
%  JACOBIAN SPARSITY PATTERN FOR SEMI-EXPLICIT DAE
%  ========================================================================

function Jpattern = compute_jacobian_pattern_semi_explicit(params)
% COMPUTE_JACOBIAN_PATTERN_SEMI_EXPLICIT Sparsity pattern for ode15s semi-explicit DAE
%
% For the semi-explicit DAE: M*y' = f(t,y) with singular mass matrix M
%
% State ordering per cell: [h, Q_l, Q_r, Q_u, Q_d, F]
%
% The Jacobian is ∂f/∂y where:
%   f_h = R - f_infil + (Q_in - Q_out) / A       (differential)
%   f_Q = Q_manning(h, h_neighbor) - Q           (algebraic constraint)
%   f_F = f_infil                                (differential)
%
% Dependencies for semi-explicit form:
%   ∂f_h/∂h:      Yes (infiltration depends on h)
%   ∂f_h/∂Q:      Yes (Q_out = Q_l + Q_r + Q_u + Q_d)
%   ∂f_h/∂F:      Yes (infiltration depends on F)
%   ∂f_h/∂Q_n:    Yes (Q_in from neighbor's Q toward me)
%
%   ∂f_Q/∂Q:      Yes (= -1, from the -Q term)
%   ∂f_Q/∂h:      Yes (Q_manning depends on self h)
%   ∂f_Q/∂h_n:    Yes (Q_manning slope depends on neighbor h)
%
%   ∂f_F/∂h:      Yes (infiltration depends on h)
%   ∂f_F/∂F:      Yes (infiltration depends on F)

    % Grid info
    grid_size = params.grid_size;
    nx = grid_size(2);
    ny = grid_size(1);
    
    if isfield(params, 'use_compressed_indexing') && params.use_compressed_indexing
        n_cells = params.valid_cells;
        grid_to_compressed = params.grid_to_compressed;
        compressed_to_grid = params.compressed_to_grid;
        use_compressed = true;
    else
        n_cells = nx * ny;
        use_compressed = false;
        grid_to_compressed = [];
        compressed_to_grid = [];
    end
    
    spc = 6;  % states per cell
    n_states = n_cells * spc;
    
    % Estimate non-zeros: ~24 per cell for semi-explicit form
    % Each cell: 6 self + up to 8 neighbor Q + up to 4 neighbor h = 18
    nnz_est = 30 * n_cells;
    I = zeros(nnz_est, 1);
    J = zeros(nnz_est, 1);
    cnt = 0;
    
    for c = 1:n_cells
        if use_compressed
            i = compressed_to_grid(c, 1);
            j = compressed_to_grid(c, 2);
        else
            i = floor((c-1)/nx) + 1;
            j = mod(c-1, nx) + 1;
        end
        
        base = (c-1) * spc;
        idx_h = base + 1;
        idx_Ql = base + 2;
        idx_Qr = base + 3;
        idx_Qu = base + 4;
        idx_Qd = base + 5;
        idx_F = base + 6;
        
        % === ∂f_h/∂* (row idx_h) ===
        % ∂f_h/∂h (infiltration depends on h)
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_h;
        % ∂f_h/∂Q_l, ∂Q_r, ∂Q_u, ∂Q_d (Q_out = sum of Q)
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Ql;
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Qr;
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Qu;
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Qd;
        % ∂f_h/∂F (infiltration depends on F)
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_F;
        
        % === ∂f_Ql/∂* (row idx_Ql) ===
        % ∂f_Ql/∂h (Q_manning_l depends on self h)
        cnt = cnt + 1; I(cnt) = idx_Ql; J(cnt) = idx_h;
        % ∂f_Ql/∂Ql (= -1 from the -Q term)
        cnt = cnt + 1; I(cnt) = idx_Ql; J(cnt) = idx_Ql;
        
        % === ∂f_Qr/∂* (row idx_Qr) ===
        cnt = cnt + 1; I(cnt) = idx_Qr; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_Qr; J(cnt) = idx_Qr;
        
        % === ∂f_Qu/∂* (row idx_Qu) ===
        cnt = cnt + 1; I(cnt) = idx_Qu; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_Qu; J(cnt) = idx_Qu;
        
        % === ∂f_Qd/∂* (row idx_Qd) ===
        cnt = cnt + 1; I(cnt) = idx_Qd; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_Qd; J(cnt) = idx_Qd;
        
        % === ∂f_F/∂* (row idx_F) ===
        cnt = cnt + 1; I(cnt) = idx_F; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_F; J(cnt) = idx_F;
        
        % === Neighbor dependencies ===
        neighbors = get_neighbor_indices(i, j, ny, nx, use_compressed, ...
            grid_to_compressed);
        
        % Left neighbor
        if neighbors(1) > 0
            nc = neighbors(1);
            base_n = (nc-1) * spc;
            % ∂f_Ql/∂h_left (slope = (wse - wse_left)/dx depends on h_left)
            cnt = cnt + 1; I(cnt) = idx_Ql; J(cnt) = base_n + 1;
            % ∂f_h/∂Qr_left (Q_in includes left neighbor's Q_r)
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 3;
        end
        
        % Right neighbor
        if neighbors(2) > 0
            nc = neighbors(2);
            base_n = (nc-1) * spc;
            % ∂f_Qr/∂h_right
            cnt = cnt + 1; I(cnt) = idx_Qr; J(cnt) = base_n + 1;
            % ∂f_h/∂Ql_right (Q_in includes right neighbor's Q_l)
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 2;
        end
        
        % Upper neighbor
        if neighbors(3) > 0
            nc = neighbors(3);
            base_n = (nc-1) * spc;
            % ∂f_Qu/∂h_up
            cnt = cnt + 1; I(cnt) = idx_Qu; J(cnt) = base_n + 1;
            % ∂f_h/∂Qd_up (Q_in includes upper neighbor's Q_d)
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 5;
        end
        
        % Lower neighbor
        if neighbors(4) > 0
            nc = neighbors(4);
            base_n = (nc-1) * spc;
            % ∂f_Qd/∂h_down
            cnt = cnt + 1; I(cnt) = idx_Qd; J(cnt) = base_n + 1;
            % ∂f_h/∂Qu_down (Q_in includes lower neighbor's Q_u)
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 4;
        end
    end
    
    % Trim and create sparse matrix
    I = I(1:cnt);
    J = J(1:cnt);
    Jpattern = sparse(I, J, ones(cnt, 1), n_states, n_states);
end


%% ========================================================================
%  JACOBIAN SPARSITY PATTERN FOR IMPLICIT DAE (kept for ode15i)
%  ========================================================================

function Jpattern = compute_jacobian_pattern_6state(params)
% COMPUTE_JACOBIAN_PATTERN_6STATE Sparsity pattern for IMPLICIT DAE (ode15i)
%
% NOTE: This function is kept for compatibility with ode15i implicit DAE solver.
%       For ode15s semi-explicit DAE, use compute_jacobian_pattern_semi_explicit().
%
% For implicit DAE: F(t, y, y') = 0
% State ordering per cell: [h, Q_l, Q_r, Q_u, Q_d, F]
%
% Dependencies:
%   dF_h/dh:      Yes (self h affects infiltration, WSE)
%   dF_h/dQ:      Yes (Q_l, Q_r, Q_u, Q_d affect mass balance)
%   dF_h/dF:      Yes (F affects infiltration)
%   dF_h/dh_n:    Yes (neighbor h affects inflow)
%   dF_h/dQ_n:    Yes (neighbor Q affects inflow)
%
%   dF_Q/dh:      Yes (h affects Manning Q)
%   dF_Q/dQ:      Yes (self Q in residual)
%   dF_Q/dh_n:    Yes (neighbor h affects slope)
%
%   dF_F/dh:      Yes (h affects infiltration)
%   dF_F/dF:      Yes (F affects infiltration)

    % Grid info
    grid_size = params.grid_size;
    nx = grid_size(2);
    ny = grid_size(1);
    
    if isfield(params, 'use_compressed_indexing') && params.use_compressed_indexing
        n_cells = params.valid_cells;
        grid_to_compressed = params.grid_to_compressed;
        compressed_to_grid = params.compressed_to_grid;
        use_compressed = true;
    else
        n_cells = nx * ny;
        use_compressed = false;
        grid_to_compressed = [];  % Not used but needed for function call
        compressed_to_grid = [];
    end
    
    spc = 6;  % states per cell
    n_states = n_cells * spc;
    
    % Estimate non-zeros: ~40 per cell (conservative)
    nnz_est = 50 * n_cells;
    I = zeros(nnz_est, 1);
    J = zeros(nnz_est, 1);
    cnt = 0;
    
    for c = 1:n_cells
        if use_compressed
            i = compressed_to_grid(c, 1);
            j = compressed_to_grid(c, 2);
        else
            i = floor((c-1)/nx) + 1;
            j = mod(c-1, nx) + 1;
        end
        
        base = (c-1) * spc;
        idx_h = base + 1;
        idx_Ql = base + 2;
        idx_Qr = base + 3;
        idx_Qu = base + 4;
        idx_Qd = base + 5;
        idx_F = base + 6;
        
        % --- df_h dependencies ---
        % df_h/dh (self)
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_h;
        % df_h/dQ_l, dQ_r, dQ_u, dQ_d (self outflows)
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Ql;
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Qr;
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Qu;
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_Qd;
        % df_h/dF (infiltration)
        cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = idx_F;
        
        % --- df_Ql dependencies ---
        % df_Ql/dh (self - Manning)
        cnt = cnt + 1; I(cnt) = idx_Ql; J(cnt) = idx_h;
        % df_Ql/dQl (self - residual)
        cnt = cnt + 1; I(cnt) = idx_Ql; J(cnt) = idx_Ql;
        
        % --- df_Qr dependencies ---
        cnt = cnt + 1; I(cnt) = idx_Qr; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_Qr; J(cnt) = idx_Qr;
        
        % --- df_Qu dependencies ---
        cnt = cnt + 1; I(cnt) = idx_Qu; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_Qu; J(cnt) = idx_Qu;
        
        % --- df_Qd dependencies ---
        cnt = cnt + 1; I(cnt) = idx_Qd; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_Qd; J(cnt) = idx_Qd;
        
        % --- df_F dependencies ---
        cnt = cnt + 1; I(cnt) = idx_F; J(cnt) = idx_h;
        cnt = cnt + 1; I(cnt) = idx_F; J(cnt) = idx_F;
        
        % --- Neighbor dependencies ---
        neighbors = get_neighbor_indices(i, j, ny, nx, use_compressed, ...
            grid_to_compressed);
        
        % Left neighbor (affects Q_l slope, h inflow via their Q_r)
        if neighbors(1) > 0
            nc = neighbors(1);
            base_n = (nc-1) * spc;
            % df_Ql/dh_left (slope depends on neighbor h)
            cnt = cnt + 1; I(cnt) = idx_Ql; J(cnt) = base_n + 1;
            % df_h/dQ_r_left (inflow from neighbor's Q_r)
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 3;
        end
        
        % Right neighbor
        if neighbors(2) > 0
            nc = neighbors(2);
            base_n = (nc-1) * spc;
            cnt = cnt + 1; I(cnt) = idx_Qr; J(cnt) = base_n + 1;
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 2;
        end
        
        % Upper neighbor
        if neighbors(3) > 0
            nc = neighbors(3);
            base_n = (nc-1) * spc;
            cnt = cnt + 1; I(cnt) = idx_Qu; J(cnt) = base_n + 1;
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 5;
        end
        
        % Lower neighbor
        if neighbors(4) > 0
            nc = neighbors(4);
            base_n = (nc-1) * spc;
            cnt = cnt + 1; I(cnt) = idx_Qd; J(cnt) = base_n + 1;
            cnt = cnt + 1; I(cnt) = idx_h; J(cnt) = base_n + 4;
        end
    end
    
    % Trim and create sparse matrix
    I = I(1:cnt);
    J = J(1:cnt);
    Jpattern = sparse(I, J, ones(cnt, 1), n_states, n_states);
end


function neighbors = get_neighbor_indices(i, j, ny, nx, use_compressed, g2c)
% GET_NEIGHBOR_INDICES Get compressed indices of 4 neighbors
%
% Returns [left, right, up, down] indices (0 if invalid/boundary)

    neighbors = zeros(4, 1);
    
    if use_compressed
        if j > 1,  neighbors(1) = g2c(i, j-1); end  % left
        if j < nx, neighbors(2) = g2c(i, j+1); end  % right
        if i > 1,  neighbors(3) = g2c(i-1, j); end  % up
        if i < ny, neighbors(4) = g2c(i+1, j); end  % down
    else
        if j > 1,  neighbors(1) = (i-1)*nx + (j-1); end
        if j < nx, neighbors(2) = (i-1)*nx + (j+1); end
        if i > 1,  neighbors(3) = (i-2)*nx + j; end
        if i < ny, neighbors(4) = i*nx + j; end
    end
end
