function pse_results = run_pse_dae(sim_params, sim_results, t_out, x_out, pse_config)
% RUN_PSE_DAE  Probabilistic State Estimation for 6-state Manning CA DAE
%
% MATHEMATICAL FRAMEWORK:
% =======================
% For BDF1 DAE: res(x^{n+1}, x^n, theta) = 0
%
% Linearized perturbation equation:
%   A * dx^{n+1} = C * dx^n + B_theta * dtheta
%
% Where:
%   A = d(res)/d(x^{n+1})  [System Jacobian]
%   C = -d(res)/d(x^n)     [Propagation matrix: diag([1,0,0,0,0,1])]
%   B_theta                [Parameter sensitivity]
%
% Covariance propagation:
%   K_xx^{n+1} = A^{-1} * K_bb * A^{-T}
%
% Where K_bb = C * K_xx^n * C^T + B_theta * K_theta * B_theta^T
%            = [temporal propagation] + [parameter injection]
%
% STATE VECTOR: x = [h, Q_l, Q_r, Q_u, Q_d, F] (6 states per cell)
%   - Differential states (h, F): Evolve via BDF, propagate variance
%   - Algebraic states (Q's): Constraints, no temporal memory
%
% DAE JACOBIAN STRUCTURE:
%   - Differential rows (h, F): I - dt*(df/dx) structure
%   - Algebraic rows (Q's): dg/dx only (NO I term)
%
% UNCERTAINTY SOURCES:
%   - Rainfall (R): CV and absolute sigma on rate [m/s]
%   - Manning n: CV, affects all Q's (creates Q-Q correlations)
%   - Infiltration (Ks, psi, eta): Log-normal Ks, affects both h and F
%
% Author: M. H. Kazma
% Date: February 2026

%% Initialize
fprintf('  Initializing FULL DAE PSE (6-state Manning CA)...\n');
fprintf('  Covariance propagation: K_xx = A^{-1} * K_bb * A^{-T}\n');

if nargin < 5 || isempty(pse_config)
	pse_config = struct();
end
if ~isfield(pse_config, 'rainfall_cv'), pse_config.rainfall_cv = 0.20; end
if ~isfield(pse_config, 'rainfall_sigma_abs'), pse_config.rainfall_sigma_abs = 0; end
if ~isfield(pse_config, 'manning_n_cv'), pse_config.manning_n_cv = 0.20; end
if ~isfield(pse_config, 'Ks_cv'), pse_config.Ks_cv = 0.20; end
if ~isfield(pse_config, 'psi_cv'), pse_config.psi_cv = 0.20; end
if ~isfield(pse_config, 'eta_cv'), pse_config.eta_cv = 0.20; end
debug_mode = isfield(pse_config, 'debug') && pse_config.debug;

% Parameter distribution configuration
param_dist = struct();
param_dist.rainfall = 'lognormal';
param_dist.manning_n = 'lognormal';
param_dist.Ks = 'lognormal';
param_dist.psi = 'lognormal';
param_dist.eta = 'truncated_gaussian';
if isfield(pse_config, 'param_distributions') && isstruct(pse_config.param_distributions)
	param_dist = merge_struct_defaults(param_dist, pse_config.param_distributions);
end

param_bounds = struct();
param_bounds.rainfall = [0, inf];
param_bounds.manning_n = [0, inf];
param_bounds.Ks = [0, inf];
param_bounds.psi = [0, inf];
param_bounds.eta = [0, 1];
if isfield(pse_config, 'param_bounds') && isstruct(pse_config.param_bounds)
	param_bounds = merge_struct_defaults(param_bounds, pse_config.param_bounds);
end

% Validate routing/formulation
if ~isfield(sim_params, 'routing_method') || ~strcmp(sim_params.routing_method, 'manning_ca')
	error('run_pse_dae supports only manning_ca routing.');
end
if ~isfield(sim_params, 'n_states_per_cell') || sim_params.n_states_per_cell < 6
	error('run_pse_dae expects 6-state formulation.');
end

% Grid info
if ~isfield(sim_params, 'grid_size')
	error('sim_params must contain grid_size.');
end
ny = sim_params.grid_size(1);
nx = sim_params.grid_size(2);

% Compressed indexing support
use_compressed = isfield(sim_params, 'use_compressed_indexing') && sim_params.use_compressed_indexing;
if use_compressed
	compressed_to_grid = sim_params.compressed_to_grid;
	grid_to_compressed = sim_params.grid_to_compressed;
	n_active = sim_params.valid_cells;
	active_mask = false(ny, nx);
	active_mask(sub2ind([ny, nx], compressed_to_grid(:,1), compressed_to_grid(:,2))) = true;
else
	n_active = nx * ny;
	active_mask = true(ny, nx);
	grid_to_compressed = [];
	compressed_to_grid = [];
end

% Parameters
Z = sim_params.Z;
dx = sim_params.dx;
dy = sim_params.dy;
cell_area = dx * dy;  % Used for neighbor inflow variance (Fix 3)

if isfield(sim_params, 'n_manning')
	n_manning = sim_params.n_manning;
	if isscalar(n_manning)
		n_manning = n_manning * ones(ny, nx);
	end
else
	error('Manning n not found in sim_params');
end

has_infiltration = isfield(sim_params, 'Ks');
if has_infiltration
	Ks = sim_params.Ks;
	eta = sim_params.eta;
	psi_f = sim_params.psi_f;
	F_min = sim_params.F_min;
end

% Outlet/boundary parameters (match compute_residual_6state defaults)
if isfield(sim_params, 'outlet_type')
	outlet_type = sim_params.outlet_type;
else
	outlet_type = 1;
end
if isfield(sim_params, 'outlet_indices')
	outlet_indices = sim_params.outlet_indices;
else
	outlet_indices = [];
end
if isfield(sim_params, 'slope_outlet')
	slope_outlet = sim_params.slope_outlet;
else
	slope_outlet = 0.01;
end
if isfield(sim_params, 'outlet_direction')
	outlet_direction = sim_params.outlet_direction;
else
	outlet_direction = 'right';
end

% Regularization / flow parameters (match compute_residual_6state defaults)
if isfield(sim_params, 'h_smooth')
	h_smooth = sim_params.h_smooth;
else
	if isfield(sim_params, 'slope_max') && sim_params.slope_max > 0.2
		h_smooth = 2e-3;
	else
		h_smooth = 1e-4;
	end
end
if isfield(sim_params, 'h_0')
	h_0 = sim_params.h_0;
else
	h_0 = 0;
end
if isfield(sim_params, 'h_flow_threshold')
	h_flow_threshold = sim_params.h_flow_threshold;
else
	h_flow_threshold = 1e-3;
end
if isfield(sim_params, 'v_limit')
	v_limit = sim_params.v_limit;
else
	v_limit = 8.0;
end

% Time discretization
nt = length(t_out);

% Optional: resample trajectory to uniform time step to reduce K_xx jumps
if isfield(pse_config, 'resample_dt_s') && pse_config.resample_dt_s > 0
	dt_resample = pse_config.resample_dt_s;
	t_uniform = (t_out(1):dt_resample:t_out(end))';
	if t_uniform(end) < t_out(end)
		t_uniform = [t_uniform; t_out(end)];
	end
	% Interpolate state trajectory onto uniform grid
	x_out = interp1(t_out(:), x_out, t_uniform, 'linear', 'extrap');
	t_out = t_uniform;
	nt = length(t_out);
end

dt_vec = diff(t_out);

% Rainfall extraction
if ~isfield(sim_results, 'R')
	error('Rainfall data not found in sim_results');
end
R_in = sim_results.R;

R_time = [];
R_pattern_scalar = 1;
R_pattern_grid = [];

use_rain_matrix_pse = false;
R_rain_matrix_grid = [];  % [ny x nx x nt] if rain_matrix mode

if isstruct(R_in)
	if isfield(R_in, 'rain_matrix')
		% Spatially distributed, time-varying rain matrix [n_cells x n_rain_steps]
		rm = R_in.rain_matrix;
		dt_rain = R_in.dt_rain;
		use_rain_matrix_pse = true;
		% Build [ny x nx x nt] grid for easy (row,col,k) lookup
		R_rain_matrix_grid = zeros(ny, nx, nt);
		for kk = 1:nt
			t_idx = min(floor(t_out(kk) / dt_rain) + 1, size(rm, 2));
			for cc = 1:n_active
				row_c = compressed_to_grid(cc, 1);
				col_c = compressed_to_grid(cc, 2);
				R_rain_matrix_grid(row_c, col_c, kk) = rm(cc, t_idx);
			end
		end
		% R_time still needed for mean calculation
		R_time = squeeze(mean(mean(R_rain_matrix_grid, 1), 2))';
		% Mean rain over active cells and time
		mask_3d = repmat(active_mask, [1, 1, nt]);
		R_vals = R_rain_matrix_grid(mask_3d);
		R_mean = mean(R_vals);
		fprintf('  PSE: rain_matrix mode (%d x %d x %d)\n', ny, nx, nt);
	elseif isfield(R_in, 'schedule')
		schedule = R_in.schedule;
		R_time = zeros(1, nt);
		for kk = 1:nt
			R_time(kk) = get_rainfall_from_schedule(schedule, t_out(kk));
		end
		if isfield(R_in, 'pattern') && ~isempty(R_in.pattern)
			[R_pattern_scalar, R_pattern_grid] = build_rain_pattern(R_in.pattern, ...
				use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
		end
	elseif isfield(R_in, 'base') && isfield(R_in, 'duration')
		base = R_in.base;
		[R_pattern_scalar, R_pattern_grid] = build_rain_pattern(base, ...
			use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
		R_time = ones(1, nt);
		R_time(t_out > R_in.duration) = 0;
	elseif isfield(R_in, 'series') && isfield(R_in, 't')
		R_time = interp1(R_in.t(:)', R_in.series(:)', t_out(:)', 'previous', 'extrap');
	else
		base = R_in;
		[R_pattern_scalar, R_pattern_grid] = build_rain_pattern(base, ...
			use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
		R_time = ones(1, nt);
	end
elseif isscalar(R_in)
	R_pattern_scalar = R_in;
	R_time = ones(1, nt);
elseif isvector(R_in) && length(R_in) == nt
	R_time = R_in(:)';
	R_pattern_scalar = 1;
else
	base = R_in;
	[R_pattern_scalar, R_pattern_grid] = build_rain_pattern(base, ...
		use_compressed, n_active, compressed_to_grid, ny, nx, active_mask);
	R_time = ones(1, nt);
end

if isempty(R_time)
	R_time = ones(1, nt);
end

if ~isempty(R_pattern_grid)
	R_mean = mean(R_pattern_grid(active_mask), 'all') * mean(R_time);
else
	R_mean = R_pattern_scalar * mean(R_time);
end

%% Extract deterministic trajectories
fprintf('  Extracting deterministic trajectories...\n');

n_states_per_cell = sim_params.n_states_per_cell;
h_traj = zeros(ny, nx, nt);
F_traj = zeros(ny, nx, nt);
Q_l_traj = zeros(ny, nx, nt);
Q_r_traj = zeros(ny, nx, nt);
Q_u_traj = zeros(ny, nx, nt);
Q_d_traj = zeros(ny, nx, nt);

for k = 1:nt
	x_k = x_out(k, :)';
	if use_compressed
		for c = 1:n_active
			i = compressed_to_grid(c,1);
			j = compressed_to_grid(c,2);
			idx_base = (c - 1) * n_states_per_cell;
			h_traj(i,j,k) = x_k(idx_base + 1);
			Q_l_traj(i,j,k) = x_k(idx_base + 2);
			Q_r_traj(i,j,k) = x_k(idx_base + 3);
			Q_u_traj(i,j,k) = x_k(idx_base + 4);
			Q_d_traj(i,j,k) = x_k(idx_base + 5);
			F_traj(i,j,k) = x_k(idx_base + 6);
		end
	else
		for idx = 1:(nx*ny)
			i = ceil(idx / nx);
			j = mod(idx - 1, nx) + 1;
			idx_base = (idx - 1) * n_states_per_cell;
			h_traj(i,j,k) = x_k(idx_base + 1);
			Q_l_traj(i,j,k) = x_k(idx_base + 2);
			Q_r_traj(i,j,k) = x_k(idx_base + 3);
			Q_u_traj(i,j,k) = x_k(idx_base + 4);
			Q_d_traj(i,j,k) = x_k(idx_base + 5);
			F_traj(i,j,k) = x_k(idx_base + 6);
		end
	end
end

%% Setup parameter variances
fprintf('  Setting up parameter variances...\n');
R_var = compute_param_variance(R_mean, pse_config.rainfall_cv, pse_config.rainfall_sigma_abs, ...
	param_dist.rainfall, param_bounds.rainfall);

cv_R_sq = pse_config.rainfall_cv^2;
sigma_R_abs_sq = pse_config.rainfall_sigma_abs^2;
cv_n_sq = pse_config.manning_n_cv^2;
cv_Ks_sq = pse_config.Ks_cv^2;
cv_psi_sq = pse_config.psi_cv^2;
cv_eta_sq = pse_config.eta_cv^2;

% Pre-compute sigma_ln_Ks ONCE (avoid isfield/sqrt/log in inner loop)
if has_infiltration
	if isfield(pse_config, 'sigma_ln_Ks')
		sig_lnKs = pse_config.sigma_ln_Ks;
	else
		sig_lnKs = sqrt(log(1 + pse_config.Ks_cv^2));
	end
else
	sig_lnKs = 0;
end

% PSE computation stride: compute every N timesteps to speed up
% Default = 1 (every timestep). User can set pse_config.stride = 5 etc.
if isfield(pse_config, 'stride') && pse_config.stride > 1
	pse_stride = pse_config.stride;
else
	pse_stride = 1;
end

% Use analytical Jacobian unless user forces numerical
use_analytical_jacobian = true;
if isfield(pse_config, 'use_numerical_jacobian') && pse_config.use_numerical_jacobian
	use_analytical_jacobian = false;
end

%% Preallocate
var_h = zeros(ny, nx, nt);
var_F = zeros(ny, nx, nt);
var_Ql = zeros(ny, nx, nt);
var_Qr = zeros(ny, nx, nt);
var_Qu = zeros(ny, nx, nt);
var_Qd = zeros(ny, nx, nt);
cov_hF = zeros(ny, nx, nt);  % Fix 1: Track h-F cross-covariance for temporal propagation

%% Step 1: FULL 6x6 DAE PSE
if use_analytical_jacobian
	fprintf('  Step 1: Computing K_xx with ANALYTICAL Jacobian (fast)...\n');
else
	fprintf('  Step 1: Computing K_xx with NUMERICAL Jacobian (slow)...\n');
end
fprintf('          Stride=%d, %d active cells, %d timesteps\n', pse_stride, n_active, nt);

% Start timing the core PSE computation (excludes preprocessing)
tic;
wall_start = tic;

active_indices = find(active_mask);
[active_rows, active_cols] = ind2sub([ny, nx], active_indices);

% Track eigenvalue cleanup usage
neg_diag_count = 0;
total_cells = 0;

if debug_mode
	min_rcond = inf;
	min_rcond_info = [NaN, NaN, NaN];
	low_rcond_count = 0;
	nonfinite_kxx_count = 0;
	max_kxx = -inf;
	max_kxx_info = [NaN, NaN, NaN, NaN];
end

% Determine which timesteps to compute PSE for
pse_timesteps = 2:pse_stride:nt;
% Always include the last timestep
if pse_timesteps(end) ~= nt
	pse_timesteps = [pse_timesteps, nt];
end
n_pse_steps = length(pse_timesteps);
pse_step_count = 0;

for ki = 1:n_pse_steps
	k = pse_timesteps(ki);
	% Previous PSE timestep index (for variance retrieval)
	if ki > 1
		k_prev_pse = pse_timesteps(ki-1);
	else
		k_prev_pse = 1;  % initial conditions
	end
	% For stride > 1, accumulate dt across skipped steps
	if pse_stride > 1 && ki > 1
		dt = sum(dt_vec(k_prev_pse:k-1));
	else
		dt = dt_vec(k-1);
	end
	
	pse_step_count = pse_step_count + 1;
	if mod(pse_step_count, 10) == 0 || pse_step_count == 1
		elapsed = toc(wall_start);
		frac = pse_step_count / n_pse_steps;
		if frac > 0.01
			eta_sec = elapsed / frac * (1 - frac);
			fprintf('    PSE step %d/%d (t=%d/%d) [%.0fs elapsed, ~%.0fs remaining]\n', ...
				pse_step_count, n_pse_steps, k, nt, elapsed, eta_sec);
		else
			fprintf('    PSE step %d/%d (t=%d/%d)\n', pse_step_count, n_pse_steps, k, nt);
		end
	end
	for a = 1:numel(active_rows)
		row = active_rows(a);
		col = active_cols(a);

		% Quick dry cell check FIRST — avoid expensive operations for dry cells
		h_curr = h_traj(row, col, k);
		dry_threshold = 1e-8;
		dry_decay_rain = 0.5;
		dry_decay_norain = 0.1;
		if isfield(pse_config, 'dry_cell_threshold'), dry_threshold = pse_config.dry_cell_threshold; end
		if isfield(pse_config, 'dry_cell_decay_rain'), dry_decay_rain = pse_config.dry_cell_decay_rain; end
		if isfield(pse_config, 'dry_cell_decay_norain'), dry_decay_norain = pse_config.dry_cell_decay_norain; end
		
		if h_curr < dry_threshold
			% --- DRY CELL FAST PATH ---
			% Propagate F variance, provide h floor; skip Jacobian entirely
			var_h_prev = var_h(row, col, k_prev_pse);
			var_F_prev = var_F(row, col, k_prev_pse);
			
			% Get rainfall for this cell
			if use_rain_matrix_pse
				R_cell = R_rain_matrix_grid(row, col, k);
			elseif ~isempty(R_pattern_grid)
				R_cell = R_pattern_grid(row, col) * R_time(k);
			else
				R_cell = R_pattern_scalar * R_time(k);
			end
			
			var_R_dry = (dt * pse_config.rainfall_cv * R_cell)^2 + ...
			            (dt * pse_config.rainfall_sigma_abs)^2;
			
			% F variance: pre-ponding, F integrates rainfall
			if has_infiltration
				Ks_cell = get_cell_param(Ks, row, col);
				var_Ks_contrib = (Ks_cell * sig_lnKs * dt)^2;
				var_F(row, col, k) = var_F_prev + var_R_dry + var_Ks_contrib;
				
				if R_cell > 0
					var_excess = var_R_dry + var_Ks_contrib;
					var_h(row, col, k) = max(var_h_prev * dry_decay_rain, var_excess);
				else
					var_h(row, col, k) = var_h_prev * dry_decay_norain;
				end
			else
				var_F(row, col, k) = var_F_prev + var_R_dry;
				if R_cell > 0
					var_h(row, col, k) = var_R_dry;
				else
					var_h(row, col, k) = var_h_prev * dry_decay_norain;
				end
			end
			
			% Q: zero for dry cells
			% var_Ql/Qr/Qu/Qd already initialized to 0
			cov_hF(row, col, k) = 0;  % Fix 1: Reset cross-covariance for dry cells
			continue;
		end

		% --- WET CELL PATH ---
		% Only extract full state for cells with water
		F_curr = F_traj(row, col, k);
		Q_l = Q_l_traj(row, col, k);
		Q_r = Q_r_traj(row, col, k);
		Q_u = Q_u_traj(row, col, k);
		Q_d = Q_d_traj(row, col, k);

		y0 = [h_curr; Q_l; Q_r; Q_u; Q_d; F_curr];

		xdot_h = (h_curr - h_traj(row, col, k_prev_pse)) / dt;
		xdot_F = (F_curr - F_traj(row, col, k_prev_pse)) / dt;

		n_cell = n_manning(row, col);
		if use_rain_matrix_pse
			R_cell = R_rain_matrix_grid(row, col, k);
		elseif ~isempty(R_pattern_grid)
			R_cell = R_pattern_grid(row, col) * R_time(k);
		else
			R_cell = R_pattern_scalar * R_time(k);
		end
		if has_infiltration
			Ks_cell = get_cell_param(Ks, row, col);
			psi_cell = get_cell_param(psi_f, row, col);
			eta_cell = get_cell_param(eta, row, col);
		else
			Ks_cell = 0; psi_cell = 0; eta_cell = 0;
		end

		% Previous variances for temporal propagation
		var_h_prev = var_h(row, col, k_prev_pse);
		var_F_prev = var_F(row, col, k_prev_pse);
		cov_hF_prev = cov_hF(row, col, k_prev_pse);  % Fix 1: Previous h-F cross-covariance

		%% Build 6x6 DAE Jacobian (System Matrix A)
		% Use ANALYTICAL Jacobian for speed
		if use_analytical_jacobian
			A_dae = build_dae_jacobian_6x6_analytical(y0, dt, ...
				n_cell, Ks_cell, psi_cell, eta_cell, ...
				dx, dy, h_0, h_smooth, h_flow_threshold, v_limit, F_min, has_infiltration);
		else
			A_dae = build_dae_jacobian_6x6_numerical(y0, row, col, k, dt, xdot_h, xdot_F, ...
				R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
				h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, ...
				Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, ...
				active_mask, F_min, has_infiltration, ...
				outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
		end

		%% Build K_bb (6x6 covariance of RHS vector b)
		%
		% MATHEMATICAL FRAMEWORK (Covariance Propagation for DAE):
		% =========================================================
		% For BDF1: res(x^{n+1}, x^n, theta) = 0
		%
		% Linearized perturbation:
		%   A * dx^{n+1} = C * dx^n + B_theta * dtheta
		%
		% Where:
		%   A = d(res)/d(x^{n+1})  [System Jacobian - built above]
		%   C = -d(res)/d(x^n)     [Propagation matrix]
		%   B_theta = -d(res)/d(theta) [Parameter sensitivity]
		%
		% For this DAE:
		%   C = diag([1, 0, 0, 0, 0, 1])
		%   - Differential states (h, F): C = 1 (x^n appears in residual)
		%   - Algebraic states (Q's): C = 0 (no Q^n term, Q determined by constraint)
		%
		% Covariance propagation:
		%   K_xx^{n+1} = A^{-1} * K_bb * A^{-T}
		%
		% Where K_bb is the covariance of the RHS:
		%   K_bb = C * K_xx^n * C^T + B_theta * K_theta * B_theta^T
		%        = [temporal propagation] + [parameter injection]
		%
		% =========================================================
		
		%% 1. PARAMETER VARIANCES
		
		% Rainfall uncertainty
		% ALWAYS include absolute floor (sigma_abs) even when R=0:
		% there is measurement/interpolation uncertainty about whether
		% rainfall is truly zero at this cell and time.
		% var(dt*R) = dt^2 * (CV*R)^2 + dt^2 * sigma_abs^2
		var_R_cell = (dt * pse_config.rainfall_cv * R_cell)^2 + (dt * pse_config.rainfall_sigma_abs)^2;
		
		% Manning n uncertainty
		% Q ~ (1/n) * h^(5/3) * sqrt(S), so dQ/dn = -Q/n
		% When Q=0 (below flow threshold), sensitivity is correctly zero:
		%   no actual flow => no uncertainty from Manning n for that face.
		% The sensitivity uses REALIZED Q (includes velocity caps, thresholds)
		% because that matches the A matrix Jacobian (same operating point).
		cv_n = pse_config.manning_n_cv;
		var_n = (cv_n * n_cell)^2;
		dQl_dn = -Q_l / n_cell;
		dQr_dn = -Q_r / n_cell;
		dQu_dn = -Q_u / n_cell;
		dQd_dn = -Q_d / n_cell;
		
		% Infiltration uncertainty (affects both h and F equations)
		% f_infil = Ks * (1 + eta*(h+psi)/F)
		if has_infiltration && F_curr > F_min
			f_infil = Ks_cell * (1 + eta_cell * (h_curr + psi_cell) / F_curr);
			
			% Log-transform for Ks (hydraulic conductivity is log-normal)
			% sig_lnKs pre-computed outside loop for speed
			df_dlnKs = f_infil;
			var_f_from_Ks = df_dlnKs^2 * sig_lnKs^2;
			
			% Standard FOSM for psi and eta
			df_dpsi = Ks_cell * eta_cell / F_curr;
			df_deta = Ks_cell * (h_curr + psi_cell) / F_curr;
			var_psi = (pse_config.psi_cv * psi_cell)^2;
			var_eta = (pse_config.eta_cv * eta_cell)^2;
			
			% Total variance of infiltration rate (not yet scaled by dt)
			var_f_rate = var_f_from_Ks + df_dpsi^2 * var_psi + df_deta^2 * var_eta;
			
			% BDF scaling: b contains dt*f, so var(b) = dt^2 * var(f)
			var_f_infil = dt^2 * var_f_rate;
		else
			var_f_infil = 0;
		end
		
		%% 2. BUILD K_bb MATRIX
		%
		% K_bb = C * K_xx^n * C^T + B_theta * K_theta * B_theta^T
		%
		% Structure of b vector in Ax = b:
		%   b_h  = h^n + dt*R - dt*f_infil  (differential)
		%   b_Ql = 0 (algebraic: res = Q - Q_manning, RHS is parameter-dependent)
		%   b_Qr = 0
		%   b_Qu = 0
		%   b_Qd = 0
		%   b_F  = F^n + dt*f_infil         (differential)
		%
		% NOTE: For algebraic equations res_Q = Q - Q_manning(h,n), 
		% the uncertainty comes from Manning n affecting Q_manning.
		% This enters through the B_theta term, not through b directly.
		% But for the K_bb formulation, we capture this as the variance
		% in the "target" value that Q must satisfy.
		
		K_bb = zeros(6, 6);
		
		% --- TEMPORAL PROPAGATION (C * K_xx^n * C^T) ---
		% C = diag([1, 0, 0, 0, 0, 1]) for [h, Ql, Qr, Qu, Qd, F]
		% Only h and F propagate from previous timestep
		% Q's are purely algebraic - no temporal memory
		K_bb(1,1) = var_h_prev;  % h^n variance propagates
		K_bb(6,6) = var_F_prev;  % F^n variance propagates
		% Fix 1: Include h-F cross-covariance from previous timestep
		% C * K_xx^n * C^T where C = diag([1,0,0,0,0,1]) includes off-diagonal
		K_bb(1,6) = cov_hF_prev;  % Temporal cross-covariance propagation
		K_bb(6,1) = cov_hF_prev;  % Symmetric
		
		% --- PARAMETER INJECTION (B_theta * K_theta * B_theta^T) ---
		
		% Rainfall affects h equation: b_h += dt*R
		K_bb(1,1) = K_bb(1,1) + var_R_cell;
		
		% Infiltration affects both h and F with OPPOSITE signs (creates correlation)
		% b_h has -dt*f (infiltration removes water from surface)
		% b_F has +dt*f (infiltration adds to cumulative)
		% Cov(b_h, b_F) = Cov(-dt*f, +dt*f) = -var(dt*f)
		K_bb(1,1) = K_bb(1,1) + var_f_infil;  % var(-dt*f) = var(dt*f)
		K_bb(6,6) = K_bb(6,6) + var_f_infil;  % var(+dt*f)
		K_bb(1,6) = K_bb(1,6) - var_f_infil;   % Add infiltration anti-correlation to temporal cross-cov
		K_bb(6,1) = K_bb(1,6);                 % Symmetric
		
		% Manning n affects ALL Q equations (creates Q-Q correlations)
		% All Q's depend on the SAME n, so they are correlated
		% Cov(Qi, Qj) = (dQi/dn) * (dQj/dn) * var(n)
		K_bb(2,2) = dQl_dn^2 * var_n;
		K_bb(3,3) = dQr_dn^2 * var_n;
		K_bb(4,4) = dQu_dn^2 * var_n;
		K_bb(5,5) = dQd_dn^2 * var_n;
		
		% Q-Q cross-covariances (same n affects all Q's)
		K_bb(2,3) = dQl_dn * dQr_dn * var_n;  K_bb(3,2) = K_bb(2,3);
		K_bb(2,4) = dQl_dn * dQu_dn * var_n;  K_bb(4,2) = K_bb(2,4);
		K_bb(2,5) = dQl_dn * dQd_dn * var_n;  K_bb(5,2) = K_bb(2,5);
		K_bb(3,4) = dQr_dn * dQu_dn * var_n;  K_bb(4,3) = K_bb(3,4);
		K_bb(3,5) = dQr_dn * dQd_dn * var_n;  K_bb(5,3) = K_bb(3,5);
		K_bb(4,5) = dQu_dn * dQd_dn * var_n;  K_bb(5,4) = K_bb(4,5);

		%% Fix 3: Neighbor inflow uncertainty
		% In the paper's global system, neighbor Q uncertainty propagates automatically.
		% For cell-local PSE, we inject variance of inflow Q from neighbors into K_bb(1,1).
		% Q_in = Q_r(left_neighbor) + Q_l(right_neighbor) + Q_d(up_neighbor) + Q_u(down_neighbor)
		% var(Q_in/A) enters the h equation via (Q_in - Q_out)/A
		var_Q_in = 0;
		if col > 1 && active_mask(row, col-1)
			var_Q_in = var_Q_in + var_Qr(row, col-1, k_prev_pse);
		end
		if col < nx && active_mask(row, col+1)
			var_Q_in = var_Q_in + var_Ql(row, col+1, k_prev_pse);
		end
		if row > 1 && active_mask(row-1, col)
			var_Q_in = var_Q_in + var_Qd(row-1, col, k_prev_pse);
		end
		if row < ny && active_mask(row+1, col)
			var_Q_in = var_Q_in + var_Qu(row+1, col, k_prev_pse);
		end
		% Scale by dt^2/cell_area^2 (inflow enters h as dt*Q_in/A in the BDF RHS)
		K_bb(1,1) = K_bb(1,1) + dt^2 * var_Q_in / cell_area^2;

		%% PSE Theorem 1: K_xx = (A^T A)^{-1} A^T K_bb A (A^T A)^{-1}
		[K_xx, rc, had_neg_diag] = compute_kxx(A_dae, K_bb);
		
		total_cells = total_cells + 1;
		if had_neg_diag
			neg_diag_count = neg_diag_count + 1;
		end

		if debug_mode
			local_max_kxx = max(abs(K_xx), [], 'all');
			if local_max_kxx > max_kxx
				max_kxx = local_max_kxx;
				max_kxx_info = [k, row, col, local_max_kxx];
			end

			if rc < min_rcond
				min_rcond = rc;
				min_rcond_info = [k, row, col];
			end
			if rc < 1e-8
				low_rcond_count = low_rcond_count + 1;
			end
			if any(~isfinite(K_xx), 'all')
				nonfinite_kxx_count = nonfinite_kxx_count + 1;
			end
		end

		% Extract variances from diagonal of K_xx
		% State order: [h, Q_l, Q_r, Q_u, Q_d, F]
		var_h(row, col, k) = max(0, K_xx(1,1));
		var_Ql(row, col, k) = max(0, K_xx(2,2));
		var_Qr(row, col, k) = max(0, K_xx(3,3));
		var_Qu(row, col, k) = max(0, K_xx(4,4));
		var_Qd(row, col, k) = max(0, K_xx(5,5));
		var_F(row, col, k) = max(0, K_xx(6,6));
		cov_hF(row, col, k) = K_xx(1,6);  % Fix 1: Store h-F cross-covariance for next timestep
	end
end

% For stride > 1, interpolate skipped timesteps
if pse_stride > 1
	stride_interp = 'linear';
	if isfield(pse_config, 'stride_interp'), stride_interp = pse_config.stride_interp; end
	fprintf('  Interpolating %d skipped timesteps (stride=%d, method=%s)...\n', ...
		nt - n_pse_steps - 1, pse_stride, stride_interp);

	if strcmp(stride_interp, 'hold')
		% Zero-order hold: copy previous computed value to skipped steps
		for ki = 2:n_pse_steps
			k_curr = pse_timesteps(ki);
			k_prev = pse_timesteps(ki-1);
			if k_curr - k_prev > 1
				for kk = (k_prev+1):(k_curr-1)
					var_h(:,:,kk) = var_h(:,:,k_prev);
					var_F(:,:,kk) = var_F(:,:,k_prev);
					var_Ql(:,:,kk) = var_Ql(:,:,k_prev);
					var_Qr(:,:,kk) = var_Qr(:,:,k_prev);
					var_Qu(:,:,kk) = var_Qu(:,:,k_prev);
					var_Qd(:,:,kk) = var_Qd(:,:,k_prev);
					cov_hF(:,:,kk) = cov_hF(:,:,k_prev);
				end
			end
		end
	elseif strcmp(stride_interp, 'pchip')
		% PCHIP interpolation per-cell (shape-preserving, C^1)
		t_pse = t_out(pse_timesteps);
		for a = 1:numel(active_rows)
			r = active_rows(a); c = active_cols(a);
			vals_h  = squeeze(var_h(r,c,pse_timesteps));
			vals_F  = squeeze(var_F(r,c,pse_timesteps));
			vals_Ql = squeeze(var_Ql(r,c,pse_timesteps));
			vals_Qr = squeeze(var_Qr(r,c,pse_timesteps));
			vals_Qu = squeeze(var_Qu(r,c,pse_timesteps));
			vals_Qd = squeeze(var_Qd(r,c,pse_timesteps));
			vals_cov = squeeze(cov_hF(r,c,pse_timesteps));
			for ki = 2:n_pse_steps
				k_curr = pse_timesteps(ki);
				k_prev = pse_timesteps(ki-1);
				if k_curr - k_prev > 1
					kk_range = (k_prev+1):(k_curr-1);
					t_interp = t_out(kk_range);
					var_h(r,c,kk_range) = max(0, interp1(t_pse, vals_h, t_interp, 'pchip'));
					var_F(r,c,kk_range) = max(0, interp1(t_pse, vals_F, t_interp, 'pchip'));
					var_Ql(r,c,kk_range) = max(0, interp1(t_pse, vals_Ql, t_interp, 'pchip'));
					var_Qr(r,c,kk_range) = max(0, interp1(t_pse, vals_Qr, t_interp, 'pchip'));
					var_Qu(r,c,kk_range) = max(0, interp1(t_pse, vals_Qu, t_interp, 'pchip'));
					var_Qd(r,c,kk_range) = max(0, interp1(t_pse, vals_Qd, t_interp, 'pchip'));
					cov_hF(r,c,kk_range) = interp1(t_pse, vals_cov, t_interp, 'pchip');
				end
			end
		end
	else
		% Linear interpolation (default)
		for ki = 2:n_pse_steps
			k_curr = pse_timesteps(ki);
			k_prev = pse_timesteps(ki-1);
			if k_curr - k_prev > 1
				for kk = (k_prev+1):(k_curr-1)
					alpha = (kk - k_prev) / (k_curr - k_prev);
					var_h(:,:,kk) = (1-alpha) * var_h(:,:,k_prev) + alpha * var_h(:,:,k_curr);
					var_F(:,:,kk) = (1-alpha) * var_F(:,:,k_prev) + alpha * var_F(:,:,k_curr);
					var_Ql(:,:,kk) = (1-alpha) * var_Ql(:,:,k_prev) + alpha * var_Ql(:,:,k_curr);
					var_Qr(:,:,kk) = (1-alpha) * var_Qr(:,:,k_prev) + alpha * var_Qr(:,:,k_curr);
					var_Qu(:,:,kk) = (1-alpha) * var_Qu(:,:,k_prev) + alpha * var_Qu(:,:,k_curr);
					var_Qd(:,:,kk) = (1-alpha) * var_Qd(:,:,k_prev) + alpha * var_Qd(:,:,k_curr);
					cov_hF(:,:,kk) = (1-alpha) * cov_hF(:,:,k_prev) + alpha * cov_hF(:,:,k_curr);  % Fix 1
				end
			end
		end
	end
end

pse_wall_time = toc(wall_start);

%% Summary and Diagnostics
fprintf('  Variance propagation complete in %.1f seconds.\n', pse_wall_time);
if use_analytical_jacobian
	fprintf('  - Method: analytical Jacobian, stride=%d\n', pse_stride);
else
	fprintf('  - Method: numerical Jacobian, stride=%d\n', pse_stride);
end
fprintf('  - h, F: Temporal propagation + parameter injection\n');
fprintf('  - Q: Parameter injection only (algebraic, no temporal memory)\n');

% PSE coverage diagnostics: report how many cells have non-zero uncertainty
last_k = nt;
var_h_final = var_h(:,:,last_k);
var_F_final = var_F(:,:,last_k);
var_Ql_final = var_Ql(:,:,last_k);
n_h_nonzero = nnz(var_h_final(active_mask) > 0);
n_F_nonzero = nnz(var_F_final(active_mask) > 0);
n_Ql_nonzero = nnz(var_Ql_final(active_mask) > 0);
fprintf('  PSE coverage at final timestep:\n');
fprintf('    var_h > 0: %d/%d cells (%.1f%%)\n', n_h_nonzero, n_active, 100*n_h_nonzero/max(n_active,1));
fprintf('    var_F > 0: %d/%d cells (%.1f%%)\n', n_F_nonzero, n_active, 100*n_F_nonzero/max(n_active,1));
fprintf('    var_Ql > 0: %d/%d cells (%.1f%%)\n', n_Ql_nonzero, n_active, 100*n_Ql_nonzero/max(n_active,1));

% Report peak variances over time
max_h_std_ever = max(sqrt(var_h(:)));
max_F_std_ever = max(sqrt(var_F(:)));
max_Ql_std_ever = max(sqrt(var_Ql(:)));
fprintf('    max h_std: %.6f m\n', max_h_std_ever);
fprintf('    max F_std: %.6f m\n', max_F_std_ever);
fprintf('    max Q_l_std: %.6f m³/s\n', max_Ql_std_ever);

%% Package results
% PSE computes ONLY covariance
% The mean/trajectory should come from the deterministic simulation (x_out)
% We do NOT store trajectory here - use x_out directly for the operating point

h_std = sqrt(var_h);
F_std = sqrt(var_F);
Q_l_std = sqrt(var_Ql);
Q_r_std = sqrt(var_Qr);
Q_u_std = sqrt(var_Qu);
Q_d_std = sqrt(var_Qd);

computation_time = toc;

pse_results = struct();

% Standard deviations (from PSE covariance diagonal)
pse_results.h_std = h_std;
pse_results.F_std = F_std;
pse_results.Q_l_std = Q_l_std;
pse_results.Q_r_std = Q_r_std;
pse_results.Q_u_std = Q_u_std;
pse_results.Q_d_std = Q_d_std;

% Full variance arrays (needed for upstream aggregation in plotting)
pse_results.var_h = var_h;
pse_results.var_F = var_F;
pse_results.var_Ql = var_Ql;
pse_results.var_Qr = var_Qr;
pse_results.var_Qu = var_Qu;
pse_results.var_Qd = var_Qd;

% Compute confidence intervals based on distribution type
if isfield(pse_config, 'distribution_type')
    dist_type = pse_config.distribution_type;
else
    dist_type = 'gaussian';  % Default
end

if isfield(pse_config, 'confidence_level')
    conf_level = pse_config.confidence_level;
else
    conf_level = 0.95;
end

fprintf('  Computing %s confidence intervals (%.0f%% level)...\n', ...
    dist_type, conf_level*100);

% Optional parameters for CI computation
ci_opts = {};
if isfield(pse_config, 'min_depth_for_gaussian')
    ci_opts = [ci_opts, {'MinDepth', pse_config.min_depth_for_gaussian}];
end
if isfield(pse_config, 'truncation_lower')
    ci_opts = [ci_opts, {'TruncLower', pse_config.truncation_lower}];
end

% Compute CIs for each state variable
[pse_results.h_lower, pse_results.h_upper] = compute_confidence_intervals(...
    h_traj, h_std, dist_type, conf_level, ci_opts{:});
pse_results.h_lower = max(pse_results.h_lower, 0);  % h >= 0 (physical)

[pse_results.F_lower, pse_results.F_upper] = compute_confidence_intervals(...
    F_traj, F_std, dist_type, conf_level, ci_opts{:});
pse_results.F_lower = max(pse_results.F_lower, 0);  % F >= 0 (physical)

% For flows, use Gaussian (can be negative, log-normal not appropriate)
[pse_results.Q_l_lower, pse_results.Q_l_upper] = compute_confidence_intervals(...
    Q_l_traj, Q_l_std, 'gaussian', conf_level);

[pse_results.Q_r_lower, pse_results.Q_r_upper] = compute_confidence_intervals(...
    Q_r_traj, Q_r_std, 'gaussian', conf_level);

[pse_results.Q_u_lower, pse_results.Q_u_upper] = compute_confidence_intervals(...
    Q_u_traj, Q_u_std, 'gaussian', conf_level);

[pse_results.Q_d_lower, pse_results.Q_d_upper] = compute_confidence_intervals(...
    Q_d_traj, Q_d_std, 'gaussian', conf_level);

% Store distribution metadata
pse_results.distribution_type = dist_type;
pse_results.confidence_level = conf_level;
pse_results.param_distributions = param_dist;
pse_results.param_bounds = param_bounds;

% Metadata
pse_results.t = t_out;
pse_results.computation_time = computation_time;
pse_results.method = sprintf('PSE Full DAE (6x6 K_xx) - %s distribution', dist_type);
fprintf('  Eigenvalue cleanup: %d/%d cells (%.1f%% had negative diagonals)\n', ...
	neg_diag_count, total_cells, 100*neg_diag_count/max(total_cells,1));
pse_results.n_cells = n_active;
pse_results.n_timesteps = nt;
pse_results.grid_size = [ny, nx];

fprintf('  PSE completed in %.2f seconds\n', computation_time);
if debug_mode
	fprintf('  Debug: min rcond(AtA)=%.2e at t_idx=%d, cell=(%d,%d)\n', ...
		min_rcond, min_rcond_info(1), min_rcond_info(2), min_rcond_info(3));
	fprintf('  Debug: rcond<1e-8 count=%d, nonfinite K_xx count=%d\n', ...
		low_rcond_count, nonfinite_kxx_count);
	fprintf('  Debug: max |K_xx|=%.3e at t_idx=%d, cell=(%d,%d)\n', ...
		max_kxx_info(4), max_kxx_info(1), max_kxx_info(2), max_kxx_info(3));
end
end

function val = get_cell_param(param, row, col)
if numel(param) > 1
	val = param(row, col);
else
	val = param;
end
end

function merged = merge_struct_defaults(base_struct, override_struct)
merged = base_struct;
fields = fieldnames(override_struct);
for ii = 1:numel(fields)
	merged.(fields{ii}) = override_struct.(fields{ii});
end
end

function var_val = compute_param_variance(mean_val, cv, sigma_abs, dist_type, bounds)
% Compute variance from mean and CV.
% Distribution choice is tracked for reporting/CI; variance uses CV definition.
var_val = (cv * mean_val)^2 + sigma_abs^2;

% Apply simple guard for bounded parameters
if strcmpi(dist_type, 'truncated_gaussian') && ~isempty(bounds)
	if isfinite(bounds(1))
		mean_val = max(mean_val, bounds(1));
	end
	if isfinite(bounds(2))
		mean_val = min(mean_val, bounds(2));
	end
	var_val = (cv * mean_val)^2 + sigma_abs^2;
end
end

function [pattern_scalar, pattern_grid] = build_rain_pattern(base, use_compressed, n_active, compressed_to_grid, ny, nx, active_mask)
pattern_grid = [];
if isscalar(base)
	pattern_scalar = base;
	return;
end

pattern_scalar = 1;
if isvector(base)
	if use_compressed && length(base) == n_active
		pattern_grid = zeros(ny, nx);
		for c = 1:n_active
			i = compressed_to_grid(c, 1);
			j = compressed_to_grid(c, 2);
			pattern_grid(i, j) = base(c);
		end
	elseif length(base) == nx*ny
		pattern_grid = reshape(base, [ny, nx]);
	else
		pattern_grid = zeros(ny, nx);
		pattern_grid(active_mask) = base(1);
	end
elseif size(base,1) == ny && size(base,2) == nx
	pattern_grid = base;
else
	pattern_grid = zeros(ny, nx);
	pattern_grid(active_mask) = base(1);
end
end

function A = build_dae_jacobian_6x6_numerical(y0, row, col, k, dt, xdot_h, xdot_F, ...
	R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
	h_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, ...
	Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, ...
	active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed)
% Build the FULL 6x6 DAE Jacobian for state x = [h, Q_l, Q_r, Q_u, Q_d, F]
%
% DAE structure:
%   Row 1 (h):    BDF residual  res_h = h^{n+1} - h^n - dt*f_h = 0
%   Row 2 (Q_l):  Algebraic     res_Ql = Q_l - Q_manning_l(h) = 0
%   Row 3 (Q_r):  Algebraic     res_Qr = Q_r - Q_manning_r(h) = 0
%   Row 4 (Q_u):  Algebraic     res_Qu = Q_u - Q_manning_u(h) = 0
%   Row 5 (Q_d):  Algebraic     res_Qd = Q_d - Q_manning_d(h) = 0
%   Row 6 (F):    BDF residual  res_F = F^{n+1} - F^n - dt*f_F = 0
%
% KEY: Differential rows have I term, algebraic rows do NOT!

cell_area = dx * dy;

% Extract states
h = max(y0(1), 0);
Q_l = y0(2);
Q_r = y0(3);
Q_u = y0(4);
Q_d = y0(5);
F = max(y0(6), F_min);

% Effective depth for Manning
h_eff = max(h - h_0, 1e-6);

% Initialize 6x6 Jacobian
A = zeros(6, 6);

%% Infiltration derivatives (exact for compute_infiltration)
if has_infiltration
	[~, df_dh, df_dF] = compute_infiltration_with_derivatives(h, F, Ks_cell, eta_cell, psi_cell, F_min);
else
	df_dh = 0;
	df_dF = 0;
end

%% Manning derivatives (numerical) to capture velocity cap + flow threshold effects
h_eps = max(1e-6, 1e-3 * max(h, 1e-6));
h_grid = h_traj(:,:,k);
h_grid(row, col) = h + h_eps;
[Qlp, Qrp, Qup, Qdp] = compute_fluxes(h_grid, row, col, n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask);

h_grid(row, col) = max(h - h_eps, 0);
[Qlm, Qrm, Qum, Qdm] = compute_fluxes(h_grid, row, col, n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask);

% Apply outlet boundary override to derivative if needed
is_outlet = false;
if ~isempty(outlet_indices)
	if use_compressed
		k_idx = grid_to_compressed(row, col);
		is_outlet = k_idx > 0 && any(outlet_indices == k_idx);
	else
		cell_idx = (row-1) * size(Z,2) + col;
		is_outlet = any(outlet_indices == cell_idx);
	end
end

if is_outlet && outlet_type == 1
	if strcmp(outlet_direction, 'right') || strcmp(outlet_direction, 'left') || ...
			strcmp(outlet_direction, 'up') || strcmp(outlet_direction, 'down')
		Qp = apply_outlet_normal(Qlp, Qrp, Qup, Qdp, h + h_eps, n_cell, dx, dy, h_smooth, h_0, h_flow_threshold, slope_outlet, outlet_direction);
		Qm = apply_outlet_normal(Qlm, Qrm, Qum, Qdm, max(h - h_eps, 0), n_cell, dx, dy, h_smooth, h_0, h_flow_threshold, slope_outlet, outlet_direction);
		Qlp = Qp(1); Qrp = Qp(2); Qup = Qp(3); Qdp = Qp(4);
		Qlm = Qm(1); Qrm = Qm(2); Qum = Qm(3); Qdm = Qm(4);
	end
end

dQ_l_dh = (Qlp - Qlm) / (2 * h_eps);
dQ_r_dh = (Qrp - Qrm) / (2 * h_eps);
dQ_u_dh = (Qup - Qum) / (2 * h_eps);
dQ_d_dh = (Qdp - Qdm) / (2 * h_eps);

dQ_out_dh = dQ_l_dh + dQ_r_dh + dQ_u_dh + dQ_d_dh;

%% Row 1: h equation (differential - BDF)
% f_h = R - f(h,F) - (Q_l + Q_r + Q_u + Q_d)/A + Q_in/A
% res_h = h^{n+1} - h^n - dt*f_h
% d(res_h)/dh = 1 - dt*(df_h/dh)

% df_h/dh = -df/dh - d(Q_out)/dh / A
df_h_dh = -df_dh - dQ_out_dh / cell_area;

% Row 1 entries (BDF structure: I - dt*(df/dx))
A(1,1) = 1 - dt * df_h_dh;           % d(res_h)/dh
A(1,2) = dt / cell_area;              % d(res_h)/dQ_l (outflow term)
A(1,3) = dt / cell_area;              % d(res_h)/dQ_r
A(1,4) = dt / cell_area;              % d(res_h)/dQ_u
A(1,5) = dt / cell_area;              % d(res_h)/dQ_d
A(1,6) = -dt * (-df_dF);              % d(res_h)/dF

%% Rows 2-5: Q equations (algebraic - NO I term!)
% Algebraic constraint: Q = Q_manning(h)
% res_Q = Q - Q_manning(h) = 0
% d(res_Q)/dQ = 1  (NOT 1 + dt*... because algebraic)
% d(res_Q)/dh = -dQ_manning/dh

% Row 2: Q_l equation
A(2,1) = -dQ_l_dh;    % d(res_Ql)/dh
A(2,2) = 1;           % d(res_Ql)/dQ_l (algebraic: just 1, no BDF scaling)

% Row 3: Q_r equation
A(3,1) = -dQ_r_dh;    % d(res_Qr)/dh
A(3,3) = 1;           % d(res_Qr)/dQ_r

% Row 4: Q_u equation
A(4,1) = -dQ_u_dh;    % d(res_Qu)/dh
A(4,4) = 1;           % d(res_Qu)/dQ_u

% Row 5: Q_d equation
A(5,1) = -dQ_d_dh;    % d(res_Qd)/dh
A(5,5) = 1;           % d(res_Qd)/dQ_d

%% Row 6: F equation (differential - BDF)
% f_F = f(h, F)  (infiltration rate)
% res_F = F^{n+1} - F^n - dt*f_F
% d(res_F)/dh = -dt*(df/dh)
% d(res_F)/dF = 1 - dt*(df/dF)

A(6,1) = -dt * df_dh;            % d(res_F)/dh
A(6,6) = 1 - dt * df_dF;         % d(res_F)/dF (BDF structure)

end

function A = build_dae_jacobian_6x6_analytical(y0, dt, ...
	n_cell, Ks_cell, psi_cell, eta_cell, ...
	dx, dy, h_0, h_smooth, h_flow_threshold, v_limit, F_min, has_infiltration)
% ANALYTICAL 6x6 DAE Jacobian —  approach
%
% Enhanced with regularization matching build_dae_jacobian_6x6.m:
%   - Mollified depth: h_reg = sqrt(h^2 + h_smooth^2)
%   - Flow factor smoothstep at h_flow_threshold
%   - Velocity cap detection (v_limit)
%
% This is ~50x faster than build_dae_jacobian_6x6 since it
% avoids compute_fluxes perturbation calls and neighbor array lookups.
%
% State order: [h, Q_l, Q_r, Q_u, Q_d, F]

cell_area = dx * dy;

h = max(y0(1), 0);
Q_l = y0(2);
Q_r = y0(3);
Q_u = y0(4);
Q_d = y0(5);
F = max(y0(6), F_min);

A = zeros(6, 6);

%% Infiltration derivatives (exact analytical)
if has_infiltration
	[~, df_dh, df_dF] = compute_infiltration_with_derivatives(h, F, Ks_cell, eta_cell, psi_cell, F_min);
else
	df_dh = 0;
	df_dF = 0;
end

%% Manning derivatives with mollified regularization
% Matches cell_fluxes physics in build_dae_jacobian_6x6.m

% Mollified effective depth
h_reg = sqrt(h^2 + h_smooth^2);
h_eff = max(h_reg - h_0, h_smooth);

% dh_eff/dh via chain rule through mollifier
if h_reg - h_0 >= h_smooth
	dh_eff_dh = h / h_reg;  % d(sqrt(h^2 + h_smooth^2))/dh
else
	dh_eff_dh = 0;  % clamped at h_smooth floor
end

% Flow factor: C^1 smoothstep transition at h_flow_threshold
if h > 0 && h_flow_threshold > 0
	t_flow = min(1.0, h / h_flow_threshold);
	flow_factor = t_flow^2 * (3 - 2*t_flow);
	if h < h_flow_threshold
		dflow_dh = 6 * t_flow * (1 - t_flow) / h_flow_threshold;
	else
		dflow_dh = 0;
	end
else
	flow_factor = 0;
	dflow_dh = 0;
end

% Velocity cap detection from Q values
% v_actual = ((Q_l+Q_r)/dy + (Q_u+Q_d)/dx) / (flow_factor * h_eff)
is_capped = false;
if flow_factor > 1e-10 && h_eff > 1e-10
	v_actual = ((Q_l + Q_r)/dy + (Q_u + Q_d)/dx) / (flow_factor * h_eff);
	is_capped = (v_actual >= v_limit * 0.99);
end

% Manning power: 5/3 for power-law, 1 when velocity-limited
if is_capped
	manning_power = 1.0;
else
	manning_power = 5/3;
end

% dQ/dh = Q * (p * dh_eff_dh / h_eff + dflow_dh / flow_factor)
if flow_factor > 1e-10
	common_factor = manning_power * dh_eff_dh / h_eff + dflow_dh / flow_factor;
else
	% Near-zero flow_factor: Q ~ 0, retain h_eff derivative only
	common_factor = manning_power * dh_eff_dh / h_eff;
end

dQ_l_dh = Q_l * common_factor;
dQ_r_dh = Q_r * common_factor;
dQ_u_dh = Q_u * common_factor;
dQ_d_dh = Q_d * common_factor;

dQ_out_dh = dQ_l_dh + dQ_r_dh + dQ_u_dh + dQ_d_dh;

%% Row 1: h equation (differential - BDF)
df_h_dh = -df_dh - dQ_out_dh / cell_area;
A(1,1) = 1 - dt * df_h_dh;
A(1,2) = dt / cell_area;
A(1,3) = dt / cell_area;
A(1,4) = dt / cell_area;
A(1,5) = dt / cell_area;
A(1,6) = dt * df_dF;  % d(res_h)/dF = -dt*(-df_dF) = dt*df_dF

%% Rows 2-5: Q algebraic constraints (NO I term)
A(2,1) = -dQ_l_dh;  A(2,2) = 1;
A(3,1) = -dQ_r_dh;  A(3,3) = 1;
A(4,1) = -dQ_u_dh;  A(4,4) = 1;
A(5,1) = -dQ_d_dh;  A(5,5) = 1;

%% Row 6: F equation (differential - BDF)
A(6,1) = -dt * df_dh;
A(6,6) = 1 - dt * df_dF;

end

function [f_h, f_F] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell, psi_cell, eta_cell, ...
	Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration)
% Differential residuals only (continuity + infiltration)
% Q's are taken from deterministic trajectory (neighbors fixed)
h = y_diff(1);
F = y_diff(2);

% Infiltration rate
if has_infiltration
	f_infil = compute_infiltration(h, F, Ks_cell, eta_cell, psi_cell, F_min);
else
	f_infil = 0;
end

% Outflow from this cell (from deterministic trajectory)
Q_out = Q_l_traj(row, col, k) + Q_r_traj(row, col, k) + Q_u_traj(row, col, k) + Q_d_traj(row, col, k);

% Inflow from neighbors (fixed from deterministic trajectory)
Q_in = 0;
ny = size(Q_l_traj, 1);
nx = size(Q_l_traj, 2);

if col > 1 && active_mask(row, col-1)
	Q_in = Q_in + Q_r_traj(row, col-1, k);
end
if col < nx && active_mask(row, col+1)
	Q_in = Q_in + Q_l_traj(row, col+1, k);
end
if row > 1 && active_mask(row-1, col)
	Q_in = Q_in + Q_d_traj(row-1, col, k);
end
if row < ny && active_mask(row+1, col)
	Q_in = Q_in + Q_u_traj(row+1, col, k);
end

cell_area = dx * dy;
f_h = xdot_h - (R_cell - f_infil + (Q_in - Q_out) / cell_area);
f_F = xdot_F - f_infil;
end

function [B_param, var_params] = local_diff_param_sens(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
	pse_config, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration)
% Parameter sensitivities for differential residuals only
[f_h0, f_F0] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell, psi_cell, eta_cell, ...
	Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);

eps_R = max(1e-10, 1e-3 * max(abs(R_cell), 1e-6));
[f_h_Rp, f_F_Rp] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell + eps_R, Ks_cell, psi_cell, eta_cell, ...
	Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
[f_h_Rm, f_F_Rm] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell - eps_R, Ks_cell, psi_cell, eta_cell, ...
	Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
dR = [(f_h_Rp - f_h_Rm)/(2*eps_R); (f_F_Rp - f_F_Rm)/(2*eps_R)];

% Manning n affects Q's but not directly the differential equations (Q's are on trajectory)
% So dn ≈ 0 for differential residuals
dn = zeros(2, 1);

if has_infiltration
	eps_Ks = max(1e-10, 1e-3 * max(abs(Ks_cell), 1e-6));
	eps_psi = max(1e-10, 1e-3 * max(abs(psi_cell), 1e-6));
	eps_eta = max(1e-10, 1e-3 * max(abs(eta_cell), 1e-6));

	[f_h_Ksp, f_F_Ksp] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell + eps_Ks, psi_cell, eta_cell, ...
		Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
	[f_h_Ksm, f_F_Ksm] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell - eps_Ks, psi_cell, eta_cell, ...
		Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
	dKs = [(f_h_Ksp - f_h_Ksm)/(2*eps_Ks); (f_F_Ksp - f_F_Ksm)/(2*eps_Ks)];

	[f_h_psip, f_F_psip] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell, psi_cell + eps_psi, eta_cell, ...
		Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
	[f_h_psim, f_F_psim] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell, psi_cell - eps_psi, eta_cell, ...
		Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
	dpsi = [(f_h_psip - f_h_psim)/(2*eps_psi); (f_F_psip - f_F_psim)/(2*eps_psi)];

	[f_h_etap, f_F_etap] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell, psi_cell, eta_cell + eps_eta, ...
		Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
	[f_h_etam, f_F_etam] = local_diff_residuals(y_diff, row, col, k, dt, xdot_h, xdot_F, R_cell, Ks_cell, psi_cell, eta_cell - eps_eta, ...
		Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, dx, dy, active_mask, F_min, has_infiltration);
	deta = [(f_h_etap - f_h_etam)/(2*eps_eta); (f_F_etap - f_F_etam)/(2*eps_eta)];

	B_param = [dR, dn, dKs, dpsi, deta];
	% Use pre-computed CV^2 factors (avoids 1.8M function calls with string comparisons)
	var_params = [
		cv_R_sq * R_cell^2 + sigma_R_abs_sq,
		cv_n_sq * n_cell^2,
		cv_Ks_sq * Ks_cell^2,
		cv_psi_sq * psi_cell^2,
		cv_eta_sq * eta_cell^2
	];
else
	B_param = [dR, dn];
	var_params = [
		cv_R_sq * R_cell^2 + sigma_R_abs_sq,
		cv_n_sq * n_cell^2
	];
end
end

function [f_infil, dfdh, dfdF, dfdKs, dfdpsi, dfdeta] = local_infiltration_derivatives(h, F, row, col, Ks, eta, psi_f, F_min)
Ks_c = get_cell_param(Ks, row, col);
eta_c = get_cell_param(eta, row, col);
psi_c = get_cell_param(psi_f, row, col);

f_infil = compute_infiltration(h, F, Ks_c, eta_c, psi_c, F_min);

% Numerical derivatives for consistency with DAE regularization
h_eps = max(1e-6, 1e-3 * max(abs(h), 1e-3));
F_eps = max(1e-6, 1e-3 * max(abs(F), 1e-3));
Ks_eps = max(1e-8, 1e-3 * max(abs(Ks_c), 1e-6));
psi_eps = max(1e-8, 1e-3 * max(abs(psi_c), 1e-6));
eta_eps = max(1e-8, 1e-3 * max(abs(eta_c), 1e-6));

f_hp = compute_infiltration(h + h_eps, F, Ks_c, eta_c, psi_c, F_min);
f_hm = compute_infiltration(max(h - h_eps, 0), F, Ks_c, eta_c, psi_c, F_min);

f_Fp = compute_infiltration(h, F + F_eps, Ks_c, eta_c, psi_c, F_min);
f_Fm = compute_infiltration(h, max(F - F_eps, 1e-8), Ks_c, eta_c, psi_c, F_min);

f_Kp = compute_infiltration(h, F, Ks_c + Ks_eps, eta_c, psi_c, F_min);
f_Km = compute_infiltration(h, F, max(Ks_c - Ks_eps, 1e-12), eta_c, psi_c, F_min);

f_psip = compute_infiltration(h, F, Ks_c, eta_c, psi_c + psi_eps, F_min);
f_psim = compute_infiltration(h, F, Ks_c, eta_c, max(psi_c - psi_eps, 0), F_min);

f_etap = compute_infiltration(h, F, Ks_c, eta_c + eta_eps, psi_c, F_min);
f_etam = compute_infiltration(h, F, Ks_c, max(eta_c - eta_eps, 1e-12), psi_c, F_min);

dfdh = (f_hp - f_hm) / (2*h_eps);
dfdF = (f_Fp - f_Fm) / (2*F_eps);
dfdKs = (f_Kp - f_Km) / (2*Ks_eps);
dfdpsi = (f_psip - f_psim) / (2*psi_eps);
dfdeta = (f_etap - f_etam) / (2*eta_eps);
end

function f_infil = compute_infiltration(h, F, Ks, eta, psi_f, F_min)
[f_infil, ~, ~] = compute_infiltration_with_derivatives(h, F, Ks, eta, psi_f, F_min);
end

function [f_infil, df_dh, df_dF] = compute_infiltration_with_derivatives(h, F, Ks, eta, psi_f, F_min)
% Exact derivatives for the same infiltration model used in compute_infiltration
h_safe = max(h, 0);
F_safe = max(F, F_min);
eta_h = eta * (h_safe + psi_f);
f_capacity = Ks * (1 + eta_h / F_safe);

h_avail = 0.005;  % 5 mm transition
x_norm = min(1.0, h_safe / h_avail);
water_factor = x_norm * x_norm * (3 - 2 * x_norm);

% d(water_factor)/dh
if h_safe <= 0
	dw_dh = 0;
elseif h_safe < h_avail
	dx_dh = 1 / h_avail;
	dw_dx = 6 * x_norm * (1 - x_norm);  % d/dx [3x^2 - 2x^3]
	dw_dh = dw_dx * dx_dh;
else
	dw_dh = 0;
end

% d(f_capacity)/dh
if h_safe <= 0
	dfc_dh = 0;
else
	dfc_dh = Ks * eta / F_safe;
end

% d(f_capacity)/dF
if F > F_min
	dfc_dF = Ks * (-eta_h / (F_safe^2));
else
	dfc_dF = 0;
end

df_dh = dfc_dh * water_factor + f_capacity * dw_dh;
df_dF = dfc_dF * water_factor;
f_infil = f_capacity * water_factor;
end

function [Q_left, Q_right, Q_up, Q_down] = compute_fluxes(h_grid, row, col, n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask)

h = h_grid(row, col);
h_safe = max(h, 0);
h_reg = sqrt(h_safe^2 + h_smooth^2);
h_eff = max(h_reg - h_0, h_smooth);

wse = Z(row, col) + h_safe;

% Neighbor WSEs: use own WSE if neighbor inactive/out of bounds
wse_left = wse; wse_right = wse; wse_up = wse; wse_down = wse;

if col > 1 && active_mask(row, col-1)
	h_left = max(h_grid(row, col-1), 0);
	wse_left = Z(row, col-1) + h_left;
end
if col < size(Z,2) && active_mask(row, col+1)
	h_right = max(h_grid(row, col+1), 0);
	wse_right = Z(row, col+1) + h_right;
end
if row > 1 && active_mask(row-1, col)
	h_up = max(h_grid(row-1, col), 0);
	wse_up = Z(row-1, col) + h_up;
end
if row < size(Z,1) && active_mask(row+1, col)
	h_down = max(h_grid(row+1, col), 0);
	wse_down = Z(row+1, col) + h_down;
end

S_left = max((wse - wse_left) / dx, 0);
S_right = max((wse - wse_right) / dx, 0);
S_up = max((wse - wse_up) / dy, 0);
S_down = max((wse - wse_down) / dy, 0);

S_total = S_left + S_right + S_up + S_down;
S_max = max([S_left, S_right, S_up, S_down]);

v_max = (1 / n_cell) * h_eff^(2/3) * sqrt(max(S_max, 0));
v_max = min(v_max, v_limit);

flow_factor = min(1.0, h_safe / h_flow_threshold);
flow_factor = flow_factor * flow_factor * (3 - 2 * flow_factor);

if S_total > 1e-12
	w_left = S_left / S_total;
	w_right = S_right / S_total;
	w_up = S_up / S_total;
	w_down = S_down / S_total;
else
	w_left = 0; w_right = 0; w_up = 0; w_down = 0;
end

Q_left = flow_factor * v_max * h_eff * dy * w_left;
Q_right = flow_factor * v_max * h_eff * dy * w_right;
Q_up = flow_factor * v_max * h_eff * dx * w_up;
Q_down = flow_factor * v_max * h_eff * dx * w_down;
end

function Q_adj = apply_outlet_normal(Q_l, Q_r, Q_u, Q_d, h_val, n_cell, dx, dy, h_smooth, h_0, h_flow_threshold, slope_outlet, outlet_direction)
% Apply outlet normal-flow override to Q components for a given h
h_safe = max(h_val, 0);
h_reg = sqrt(h_safe^2 + h_smooth^2);
h_eff = max(h_reg - h_0, h_smooth);

flow_factor = min(1.0, h_safe / h_flow_threshold);
flow_factor = flow_factor * flow_factor * (3 - 2 * flow_factor);

if strcmp(outlet_direction, 'right')
	Q_normal = flow_factor * (dy * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
	Q_r = Q_normal;
elseif strcmp(outlet_direction, 'down')
	Q_normal = flow_factor * (dx * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
	Q_d = Q_normal;
elseif strcmp(outlet_direction, 'left')
	Q_normal = flow_factor * (dy * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
	Q_l = Q_normal;
elseif strcmp(outlet_direction, 'up')
	Q_normal = flow_factor * (dx * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
	Q_u = Q_normal;
end

Q_adj = [Q_l, Q_r, Q_u, Q_d];
end

function Jx = local_jacobian(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
	h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed)
% ANALYTICAL Jacobian of local residuals (6x6)
% Matches run_manning_ca_bdf_newton.m approach using dQ/dh = (5/3)*Q/h_eff

Jx = zeros(6, 6);
cell_area = dx * dy;

% Extract states
h = max(y0(1), 0);
Q_l = y0(2);
Q_r = y0(3);
Q_u = y0(4);
Q_d = y0(5);
F = max(y0(6), F_min);

%% Infiltration derivatives (exact for compute_infiltration)
if has_infiltration
	[~, df_dh, df_dF] = compute_infiltration_with_derivatives(h, F, Ks_cell, eta_cell, psi_cell, F_min);
else
	df_dh = 0;
	df_dF = 0;
end

%% Analytical Manning derivatives: dQ/dh = (5/3)*Q/h_eff
h_eff = max(h - h_0, 0);
if h_eff > 1e-6
	dQ_l_dh = (5/3) * Q_l / h_eff;
	dQ_r_dh = (5/3) * Q_r / h_eff;
	dQ_u_dh = (5/3) * Q_u / h_eff;
	dQ_d_dh = (5/3) * Q_d / h_eff;
else
	dQ_l_dh = 0;
	dQ_r_dh = 0;
	dQ_u_dh = 0;
	dQ_d_dh = 0;
end

dQ_out_dh = dQ_l_dh + dQ_r_dh + dQ_u_dh + dQ_d_dh;

%% Row 1: ∂res_h/∂[h, Q_l, Q_r, Q_u, Q_d, F]
% res_h = h - h_old - dt*(R - f + (Q_in - Q_out)/A)
% ∂res_h/∂h = 1 + dt*df/dh + dt*dQ_out/dh/A
Jx(1,1) = 1 + dt * df_dh + dt * dQ_out_dh / cell_area;

% ∂res_h/∂Q_i = +dt/A (positive because Q_out appears with minus sign)
Jx(1,2) = dt / cell_area;  % Q_l
Jx(1,3) = dt / cell_area;  % Q_r
Jx(1,4) = dt / cell_area;  % Q_u
Jx(1,5) = dt / cell_area;  % Q_d

% ∂res_h/∂F = dt*df/dF
Jx(1,6) = dt * df_dF;

%% Rows 2-5: ∂res_Q/∂[h, Q_l, Q_r, Q_u, Q_d, F]
% res_Q_i = Q_i - Q_manning_i(h)
% ∂res_Q_i/∂Q_i = 1
Jx(2,2) = 1;  % ∂res_Q_l/∂Q_l
Jx(3,3) = 1;  % ∂res_Q_r/∂Q_r
Jx(4,4) = 1;  % ∂res_Q_u/∂Q_u
Jx(5,5) = 1;  % ∂res_Q_d/∂Q_d

% ∂res_Q_i/∂h = -dQ_manning_i/dh
Jx(2,1) = -dQ_l_dh;  % ∂res_Q_l/∂h
Jx(3,1) = -dQ_r_dh;  % ∂res_Q_r/∂h
Jx(4,1) = -dQ_u_dh;  % ∂res_Q_u/∂h
Jx(5,1) = -dQ_d_dh;  % ∂res_Q_d/∂h

%% Row 6: ∂res_F/∂[h, Q_l, Q_r, Q_u, Q_d, F]
% res_F = F - F_old - dt*f
% ∂res_F/∂h = -dt*df/dh
Jx(6,1) = -dt * df_dh;

% ∂res_F/∂F = 1 - dt*df/dF
Jx(6,6) = 1 - dt * df_dF;

end

function [B_param, var_params] = local_param_sens(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
	pse_config, h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed)
% Sensitivities of residuals w.r.t uncertain parameters
f0 = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
	h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);

eps_R = max(1e-10, 1e-3 * max(abs(R_cell), 1e-6));
fR_p = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell + eps_R, n_cell, Ks_cell, psi_cell, eta_cell, ...
	h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
fR_m = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell - eps_R, n_cell, Ks_cell, psi_cell, eta_cell, ...
	h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
dR = (fR_p - fR_m) / (2 * eps_R);

eps_n = max(1e-8, 1e-3 * max(abs(n_cell), 1e-3));
fn_p = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell + eps_n, Ks_cell, psi_cell, eta_cell, ...
	h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
fn_m = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell - eps_n, Ks_cell, psi_cell, eta_cell, ...
	h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
dn = (fn_p - fn_m) / (2 * eps_n);

	if has_infiltration
		eps_Ks = max(1e-10, 1e-3 * max(abs(Ks_cell), 1e-6));
		eps_psi = max(1e-10, 1e-3 * max(abs(psi_cell), 1e-6));
		eps_eta = max(1e-10, 1e-3 * max(abs(eta_cell), 1e-6));

		fKs_p = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell + eps_Ks, psi_cell, eta_cell, ...
			h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
			outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
		fKs_m = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell - eps_Ks, psi_cell, eta_cell, ...
			h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
			outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
		dKs = (fKs_p - fKs_m) / (2 * eps_Ks);

		fpsi_p = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell + eps_psi, eta_cell, ...
			h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
			outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
		fpsi_m = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell - eps_psi, eta_cell, ...
			h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
			outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
		dpsi = (fpsi_p - fpsi_m) / (2 * eps_psi);

		feta_p = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell, eta_cell + eps_eta, ...
			h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
			outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
		feta_m = local_residuals(y0, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell, eta_cell - eps_eta, ...
			h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
			outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed);
		deta = (feta_p - feta_m) / (2 * eps_eta);

		B_param = [dR, dn, dKs, dpsi, deta];
		% Use pre-computed CV^2 factors (avoids function calls with string comparisons)
		var_params = [
			cv_R_sq * R_cell^2 + sigma_R_abs_sq,
			cv_n_sq * n_cell^2,
			cv_Ks_sq * Ks_cell^2,
			cv_psi_sq * psi_cell^2,
			cv_eta_sq * eta_cell^2
		];
	else
	B_param = [dR, dn];
	var_params = [
			cv_R_sq * R_cell^2 + sigma_R_abs_sq,
			cv_n_sq * n_cell^2
	];
end
end

function f = local_residuals(y, row, col, k, dt, xdot_h, xdot_F, R_cell, n_cell, Ks_cell, psi_cell, eta_cell, ...
	h_traj, F_traj, Q_l_traj, Q_r_traj, Q_u_traj, Q_d_traj, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask, F_min, has_infiltration, ...
	outlet_type, outlet_indices, outlet_direction, slope_outlet, use_compressed, grid_to_compressed)
% Local residuals for one cell (6 equations)
h = y(1);
Q_l = y(2);
Q_r = y(3);
Q_u = y(4);
Q_d = y(5);
F = y(6);

h_grid = h_traj(:,:,k);
h_grid(row, col) = h;

% Inflow from neighbors (fixed from deterministic Q)
Q_in = 0;
if col > 1 && active_mask(row, col-1)
	Q_in = Q_in + Q_r_traj(row, col-1, k);
end
if col < size(h_grid,2) && active_mask(row, col+1)
	Q_in = Q_in + Q_l_traj(row, col+1, k);
end
if row > 1 && active_mask(row-1, col)
	Q_in = Q_in + Q_d_traj(row-1, col, k);
end
if row < size(h_grid,1) && active_mask(row+1, col)
	Q_in = Q_in + Q_u_traj(row+1, col, k);
end

Q_out = Q_l + Q_r + Q_u + Q_d;

if has_infiltration
	f_infil = compute_infiltration(h, F, Ks_cell, eta_cell, psi_cell, F_min);
else
	f_infil = 0;
end

cell_area = dx * dy;
f1 = xdot_h - (R_cell - f_infil + (Q_in - Q_out) / cell_area);
f6 = xdot_F - f_infil;

[Qlc, Qrc, Quc, Qdc] = compute_fluxes(h_grid, row, col, n_cell, Z, dx, dy, h_smooth, h_0, h_flow_threshold, v_limit, active_mask);

f2 = Q_l - Qlc;
f3 = Q_r - Qrc;
f4 = Q_u - Quc;
f5 = Q_d - Qdc;

% Apply outlet boundary conditions (match compute_residual_6state)
is_outlet = false;
if ~isempty(outlet_indices)
	if use_compressed
		k_idx = grid_to_compressed(row, col);
		is_outlet = k_idx > 0 && any(outlet_indices == k_idx);
	else
		cell_idx = (row-1) * size(Z,2) + col;
		is_outlet = any(outlet_indices == cell_idx);
	end
end

if is_outlet && outlet_type == 1
	h_safe = max(h, 0);
	h_reg = sqrt(h_safe^2 + h_smooth^2);
	h_eff = max(h_reg - h_0, h_smooth);
	flow_factor = min(1.0, h_safe / h_flow_threshold);
	flow_factor = flow_factor * flow_factor * (3 - 2 * flow_factor);

	if strcmp(outlet_direction, 'right')
		Q_normal = flow_factor * (dy * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
		f3 = Q_r - Q_normal;
	elseif strcmp(outlet_direction, 'down')
		Q_normal = flow_factor * (dx * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
		f5 = Q_d - Q_normal;
	elseif strcmp(outlet_direction, 'left')
		Q_normal = flow_factor * (dy * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
		f2 = Q_l - Q_normal;
	elseif strcmp(outlet_direction, 'up')
		Q_normal = flow_factor * (dx * h_eff / n_cell) * h_eff^(2/3) * sqrt(slope_outlet);
		f4 = Q_u - Q_normal;
	end
end

f = [f1; f2; f3; f4; f5; f6];
end

function [Kxx, rc, had_neg_diag] = compute_kxx(A, Kbb)
% Covariance propagation through linear system Ax = b
%
% Given: Ax = b where b has covariance K_bb
% Then:  x = A^{-1} b
% And:   K_xx = A^{-1} K_bb A^{-T}
%
% For non-square or rank-deficient A, use pseudo-inverse form:
%   K_xx = (A^T A)^{-1} A^T K_bb A (A^T A)^{-1}
%
% REGULARIZATION: Add small diagonal term to handle rank deficiency
%   when cells are dry (h -> 0) or constraints degenerate

% Ensure K_bb is symmetric
Kbb = 0.5 * (Kbb + Kbb');

% Check condition number of A directly (more relevant than A^T A)
rc = rcond(A);

if rc > 1e-12 && ~isnan(rc)
    % Well-conditioned square matrix: use direct formula
    % K_xx = A^{-1} K_bb A^{-T}
    A_inv = A \ eye(size(A));
    Kxx = A_inv * Kbb * A_inv';
else
    % Ill-conditioned or singular: use regularized pseudo-inverse
    AtA = A' * A;
    
    % Tikhonov regularization
    lambda = 1e-10 * max(max(abs(diag(AtA))), 1e-16);
    AtA_reg = AtA + lambda * eye(size(AtA));
    AtA_inv = AtA_reg \ eye(size(AtA_reg));
    
    AtA_inv_At = AtA_inv * A';
    Kxx = AtA_inv_At * Kbb * AtA_inv_At';
end

% Ensure output is symmetric
Kxx = 0.5 * (Kxx + Kxx');

% Only do expensive eigenvalue cleanup if we detect negative diagonal elements
% (which would indicate potential negative eigenvalues)
had_neg_diag = any(diag(Kxx) < 0);
if had_neg_diag
    % Clamp small negative eigenvalues to zero (numerical cleanup)
    [V, D] = eig(Kxx);
    D = diag(max(diag(D), 0));
    Kxx = V * D * V';
    Kxx = 0.5 * (Kxx + Kxx');  % Re-symmetrize after eigen decomposition
end
end
