function [t_out, x_out, results] = run_manning_ca_bdf_newton(params, R_input, duration, varargin)
%
% Outputs full state per cell: [h, Q_l, Q_r, Q_u, Q_d, F]

%% Parse inputs
p = inputParser;
addParameter(p, 'dt', 1.0);
addParameter(p, 'BDFOrder', 1);
addParameter(p, 'NewtonTol', 1e-6);      % RMS tolerance
addParameter(p, 'MaxNewtonIter', 15);
addParameter(p, 'MaxLineSearch', 10);
addParameter(p, 'DisplayProgress', true);
addParameter(p, 'OutputInterval', 60);
addParameter(p, 'MaxStep', []);          % not used (kept for compatibility)
parse(p, varargin{:});

dt = p.Results.dt;
bdf_order = p.Results.BDFOrder;
newton_tol = p.Results.NewtonTol;
max_newton_iter = p.Results.MaxNewtonIter;
max_ls = p.Results.MaxLineSearch;
display_progress = p.Results.DisplayProgress;
output_interval = p.Results.OutputInterval;

%% Extract params
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

% Rainfall vector
if isscalar(R_input)
    R = R_input * ones(n_cells, 1);
else
    R = R_input(:);
    if numel(R) ~= n_cells
        error('R_input length mismatch: expected %d cells', n_cells);
    end
end

% Infiltration
if isfield(params, 'Ks')
    Ks = params.Ks;
    eta = params.eta;
    psi_f = params.psi_f;
    F_min = params.F_min;
    F_initial = F_min;
    if isfield(params, 'F_initial')
        F_initial = max(params.F_initial, F_min);
    end
    use_infiltration = true;
else
    use_infiltration = false;
    Ks = [];
    eta = [];
    psi_f = [];
end

% Outlet
if isfield(params, 'outlet_direction'), outlet_direction = params.outlet_direction; else, outlet_direction = 'down'; end
if isfield(params, 'slope_outlet'), slope_outlet = params.slope_outlet; else, slope_outlet = 0.01; end
if isfield(params, 'outlet_indices'), outlet_indices = params.outlet_indices; else, outlet_indices = []; end
is_outlet_cell = false(n_cells, 1);
is_outlet_cell(outlet_indices) = true;

%% BDF coefficients
if bdf_order == 1
    alpha = [1]; beta = 1;       % backward Euler
elseif bdf_order == 2
    alpha = [4/3, -1/3]; beta = 2/3;
else
    error('Only BDF order 1 and 2 supported');
end

%% Reduced state: 2 per cell [h, F]
spc = 2;
n_states_red = n_cells * spc;

% Full output state: 6 per cell [h, Ql, Qr, Qu, Qd, F]
spc_full = 6;
n_states_full = n_cells * spc_full;

%% Initialize reduced state
xHF = zeros(n_states_red, 1);
for c = 1:n_cells
    b2 = (c-1)*2;
    xHF(b2 + 1) = 0;            % h
    xHF(b2 + 2) = F_initial;    % F (antecedent soil moisture)
end
xHF_prev = xHF;

%% Outputs
n_steps = ceil(duration / dt);
n_outputs = floor(n_steps / output_interval) + 1;
t_out = zeros(n_outputs, 1);
x_out = zeros(n_outputs, n_states_full);

t = 0;
output_idx = 1;

[t_out(output_idx), x_out(output_idx,:)] = pack_full_state(t, xHF, params, ...
    use_compressed, n_cells, grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, Z, n_manning, h_0, ...
    outlet_direction, slope_outlet, is_outlet_cell);
output_idx = output_idx + 1;

if display_progress
    fprintf('dt=%.3g, BDF%d, NewtonTol(RMS)=%.1e, MaxNewton=%d\n', dt, bdf_order, newton_tol, max_newton_iter);
end

tic;
total_newton = 0;
total_outflow = 0;  % Track cumulative boundary outflow

%% Time-stepping
for step = 1:n_steps
    t_new = t + dt;
    dt_actual = t_new - t;

    xHF_new = xHF;  % initial guess

    % Newton iterations
    for it = 1:max_newton_iter
        % residual + Jacobian
        [res, J] = residual_and_jacobian_HF(xHF_new, xHF, xHF_prev, step, dt, bdf_order, alpha, beta, ...
            R, params, use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
            nx, ny, dx, dy, cell_area, Z, n_manning, h_0, F_min, ...
            Ks, eta, psi_f, use_infiltration, outlet_direction, slope_outlet, is_outlet_cell);

        res_rms = norm(res) / sqrt(length(res));

        if res_rms < newton_tol
            break;
        end

        % Newton step
        dxHF = -(J \ res);

        % Increment RMS check (scale with state size)
        dx_rms = norm(dxHF) / sqrt(length(dxHF));
        if dx_rms < newton_tol
            xHF_new = xHF_new + dxHF;
            break;
        end

        % Line search damping (Armijo-like)
        lambda = 1.0;
        res0 = res_rms;
        x_try = xHF_new;
        ok = false;

        for ls = 1:max_ls
            x_try = xHF_new + lambda * dxHF;
            res_try = residual_only_HF(x_try, xHF, xHF_prev, step, dt, bdf_order, alpha, beta, ...
                R, params, use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
                nx, ny, dx, dy, cell_area, Z, n_manning, h_0, F_min, ...
                Ks, eta, psi_f, use_infiltration, outlet_direction, slope_outlet, is_outlet_cell);

            res_try_rms = norm(res_try) / sqrt(length(res_try));

            if res_try_rms < (1 - 1e-4*lambda) * res0
                ok = true;
                break;
            end
            lambda = 0.5 * lambda;
        end

        if ~ok
            % fallback: accept damped step anyway (tiny) or break
            xHF_new = xHF_new + (0.1)*dxHF;
        else
            xHF_new = x_try;
        end
    end

    total_newton = total_newton + it;

    % Post-step physical bounds (safe)
    for c = 1:n_cells
        b2 = (c-1)*2;
        xHF_new(b2+1) = max(xHF_new(b2+1), 0);      % h >= 0
        xHF_new(b2+2) = max(xHF_new(b2+2), F_min);  % F >= F_min
    end

    % Compute boundary outflow from converged state
    Q_boundary = compute_boundary_outflow_newton(xHF_new, params, use_compressed, n_cells, ...
        grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, Z, n_manning, h_0, ...
        outlet_direction, slope_outlet, is_outlet_cell);
    
    % Track total outflow
    total_outflow = total_outflow + Q_boundary * dt_actual;
    
    % shift history
    if bdf_order == 2
        xHF_prev = xHF;
    end
    xHF = xHF_new;
    t = t_new;

    % output
    if mod(step, output_interval) == 0 || t >= duration
        [t_out(output_idx), x_out(output_idx,:)] = pack_full_state(t, xHF, params, ...
            use_compressed, n_cells, grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, Z, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell);
        if display_progress
            h_vec = xHF(1:2:end);
            F_vec = xHF(2:2:end);
            fprintf('t=%7.1fs (%.1fmin): h_max=%.4fm, h_mean=%.3fmm, F_max=%.1fmm, Newton=%d\n', ...
                t, t/60, max(h_vec), mean(h_vec)*1000, max(F_vec)*1000, it);
        end
        output_idx = output_idx + 1;
        if output_idx > n_outputs, break; end
    end
end

comp_time = toc;

% Trim output arrays
t_out = t_out(1:output_idx-1);
x_out = x_out(1:output_idx-1,:);

%% Water balance
h_final = xHF(1:2:end);
F_final = xHF(2:2:end);

total_rainfall = mean(R) * duration * n_cells * cell_area;
total_ponding  = sum(h_final) * cell_area;
total_infil    = sum(F_final - F_initial) * cell_area;

% Use tracked boundary outflow for runoff
total_runoff   = total_outflow;

% Water balance check: Input = Outputs + Storage
% Rainfall = Infiltration + Runoff + Ponding + Error
balance_error  = total_rainfall - (total_infil + total_runoff + total_ponding);

%% Store results
results.h_final = h_final;
results.F_final = F_final;
results.R = R;  % Store rainfall intensity for PSE
results.computation_time = comp_time;
results.total_newton_iter = total_newton;
results.avg_newton_per_step = total_newton / max(1, n_steps);

results.water_balance.rainfall = total_rainfall;
results.water_balance.ponding = total_ponding;
results.water_balance.infiltration = total_infil;
results.water_balance.runoff = total_runoff;
results.water_balance.error = balance_error;
results.water_balance.error_pct = 100 * abs(balance_error) / max(total_rainfall, 1);

if display_progress
    fprintf('========================================\n');
    fprintf('ROUTE 1 done. Time: %.2fs, Avg Newton/step: %.2f\n', comp_time, results.avg_newton_per_step);
    fprintf('\nWater Balance:\n');
    fprintf('  Rainfall:     %10.2f m³\n', total_rainfall);
    fprintf('  Infiltration: %10.2f m³ (%.1f%%)\n', total_infil, 100*total_infil/total_rainfall);
    fprintf('  Ponding:      %10.2f m³ (%.1f%%)\n', total_ponding, 100*total_ponding/total_rainfall);
    fprintf('  Runoff:       %10.2f m³ (%.1f%%)\n', total_runoff, 100*total_runoff/total_rainfall);
    fprintf('  Error:        %10.2f m³ (%.4f%%)\n', balance_error, results.water_balance.error_pct);
    fprintf('\nFinal State:\n');
    fprintf('  Water depth:  min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
        min(h_final(h_final>0))*1000, max(h_final)*1000, mean(h_final)*1000);
    fprintf('  Infiltration: min=%.2f mm, max=%.2f mm, mean=%.2f mm\n', ...
        min(F_final)*1000, max(F_final)*1000, mean(F_final)*1000);
end

end

%% ========================================================================
%  Residual
%  ========================================================================
function [res, J] = residual_and_jacobian_HF(xHF_new, xHF, xHF_prev, step, dt, bdf_order, alpha, beta, ...
    R, params, use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
    nx, ny, dx, dy, cell_area, Z, n_manning, h_0, F_min, ...
    Ks, eta, psi_f, use_infiltration, outlet_direction, slope_outlet, is_outlet_cell)

n = n_cells*2;
res = zeros(n,1);

% Triplets for sparse J
% rough nnz: per cell about (row h: up to 6 entries, row F: 2 entries) => ~8*n_cells
I = zeros(8*n_cells,1);
Jc = zeros(8*n_cells,1);
V = zeros(8*n_cells,1);
nz = 0;

heps = 1e-9;   % smooth regularization for h_eff in powers/denoms

for c = 1:n_cells
    % map to grid
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1)/nx) + 1;
        j = mod(c-1, nx) + 1;
    end

    b2 = (c-1)*2;
    h  = xHF_new(b2+1);
    F  = xHF_new(b2+2);

    h_old = xHF(b2+1);
    F_old = xHF(b2+2);

    if bdf_order==2 && step>1
        h_prev = xHF_prev(b2+1);
        F_prev = xHF_prev(b2+2);
    else
        h_prev = 0; F_prev = 0;
    end

    % Neighbor h (for coupling) + indices
    [hL, idxL] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'left',  c);
    [hR, idxR] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'right', c);
    [hU, idxU] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'up',    c);
    [hD, idxD] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'down',  c);

    % Infiltration f(h,F) and derivatives
    [f_infil, df_dh, df_dF] = infiltration_smooth(h, F, dt, i, j, Ks, eta, psi_f, F_min, use_infiltration);

    % Manning flows Q(h, neighbor h) and approximate derivatives wrt own and neighbor h
    is_outlet = is_outlet_cell(c);
    nloc = get_manning(n_manning, i, j);

    [Q, dQdh_own, dQdh_nb] = manning_Q_and_dQ(h, hL, hR, hU, hD, Z, i, j, nx, ny, dx, dy, nloc, h_0, ...
        outlet_direction, slope_outlet, is_outlet, heps);

    Ql = Q(1); Qr = Q(2); Qu = Q(3); Qd = Q(4);

    % inflows from neighbors use THEIR outflow toward us; with elimination we model net via continuity:
    % Here we use the same convention as your original:
    % Q_out = Ql+Qr+Qu+Qd (out of cell)
    % Q_in  = neighbor flows into cell: left neighbor's Qr, right neighbor's Ql, up neighbor's Qd, down neighbor's Qu
    % We must compute those using neighbor cell Manning too to remain consistent.
    % For speed we approximate by using symmetry: inflow depends on neighbor h and our h via that neighbor's manning;
    % We compute neighbor -> us flux by calling manning for neighbor cell (4 calls).
    Q_in = 0;
    dQin_dh_own = 0;     % derivative of Qin wrt our h
    dQin_dh_nb = struct('L',0,'R',0,'U',0,'D',0); % wrt neighbor h

    % Left neighbor contributes its Qr into us
    if idxL > 0
        [Qn, dOwn, dNb] = manning_Q_neighbor_into_us(xHF_new, idxL, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromLeft');
        Q_in = Q_in + Qn;
        dQin_dh_nb.L = dQin_dh_nb.L + dOwn;   % neighbor's own h affects its outgoing to us
        dQin_dh_own  = dQin_dh_own  + dNb;    % our h affects neighbor slope term (approx)
    end

    % Right neighbor contributes its Ql into us
    if idxR > 0
        [Qn, dOwn, dNb] = manning_Q_neighbor_into_us(xHF_new, idxR, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromRight');
        Q_in = Q_in + Qn;
        dQin_dh_nb.R = dQin_dh_nb.R + dOwn;
        dQin_dh_own  = dQin_dh_own  + dNb;
    end

    % Up neighbor contributes its Qd into us
    if idxU > 0
        [Qn, dOwn, dNb] = manning_Q_neighbor_into_us(xHF_new, idxU, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromUp');
        Q_in = Q_in + Qn;
        dQin_dh_nb.U = dQin_dh_nb.U + dOwn;
        dQin_dh_own  = dQin_dh_own  + dNb;
    end

    % Down neighbor contributes its Qu into us
    if idxD > 0
        [Qn, dOwn, dNb] = manning_Q_neighbor_into_us(xHF_new, idxD, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromDown');
        Q_in = Q_in + Qn;
        dQin_dh_nb.D = dQin_dh_nb.D + dOwn;
        dQin_dh_own  = dQin_dh_own  + dNb;
    end

    Q_out = Ql + Qr + Qu + Qd;
    Q_net = Q_in - Q_out;

    % dh/dt
    dh_dt = R(c) - f_infil + Q_net / cell_area;

    if bdf_order == 1
        res_h = h - h_old - dt*beta*dh_dt;
        res_F = F - F_old - dt*beta*f_infil;
    else
        res_h = h - alpha(1)*h_old - alpha(2)*h_prev - dt*beta*dh_dt;
        res_F = F - alpha(1)*F_old - alpha(2)*F_prev - dt*beta*f_infil;
    end

    res(b2+1) = res_h;
    res(b2+2) = res_F;

    % Jacobian entries
    % res_h = h - (...) - dt*beta*(R - f + (Qin-Qout)/A)
    % drh/dh = 1 - dt*beta*( -df/dh + (dQin/dh - dQout/dh)/A )
    %        = 1 + dt*beta*df/dh - dt*beta*(dQin_dh - dQout_dh)/A
    dQout_dh_own = sum(dQdh_own);
    dQout_dh_L = dQdh_nb(1);
    dQout_dh_R = dQdh_nb(2);
    dQout_dh_U = dQdh_nb(3);
    dQout_dh_D = dQdh_nb(4);

    drh_dh = 1 + dt*beta*df_dh - dt*beta*( (dQin_dh_own - dQout_dh_own)/cell_area );
    drh_dF = dt*beta*df_dF;  % because -dt*beta*(R - f + ...) gives +dt*beta*df/dF

    % Fill row for res_h
    nz=nz+1; I(nz)=b2+1; Jc(nz)=b2+1; V(nz)=drh_dh;
    nz=nz+1; I(nz)=b2+1; Jc(nz)=b2+2; V(nz)=drh_dF;

    % Neighbor h couplings in res_h (through Qin and Qout)
    if idxL>0
        dQnet_dhL = dQin_dh_nb.L - dQout_dh_L;
        nz=nz+1; I(nz)=b2+1; Jc(nz)=(idxL-1)*2+1; V(nz)= -dt*beta*(dQnet_dhL/cell_area);
    end
    if idxR>0
        dQnet_dhR = dQin_dh_nb.R - dQout_dh_R;
        nz=nz+1; I(nz)=b2+1; Jc(nz)=(idxR-1)*2+1; V(nz)= -dt*beta*(dQnet_dhR/cell_area);
    end
    if idxU>0
        dQnet_dhU = dQin_dh_nb.U - dQout_dh_U;
        nz=nz+1; I(nz)=b2+1; Jc(nz)=(idxU-1)*2+1; V(nz)= -dt*beta*(dQnet_dhU/cell_area);
    end
    if idxD>0
        dQnet_dhD = dQin_dh_nb.D - dQout_dh_D;
        nz=nz+1; I(nz)=b2+1; Jc(nz)=(idxD-1)*2+1; V(nz)= -dt*beta*(dQnet_dhD/cell_area);
    end

    % res_F = F - (...) - dt*beta*f
    drF_dh = -dt*beta*df_dh;
    drF_dF =  1 -dt*beta*df_dF;

    nz=nz+1; I(nz)=b2+2; Jc(nz)=b2+1; V(nz)=drF_dh;
    nz=nz+1; I(nz)=b2+2; Jc(nz)=b2+2; V(nz)=drF_dF;
end

I = I(1:nz); Jc = Jc(1:nz); V = V(1:nz);
J = sparse(I,Jc,V,n,n);

end

%% Residual only (for line search)
function res = residual_only_HF(xHF_new, xHF, xHF_prev, step, dt, bdf_order, alpha, beta, ...
    R, params, use_compressed, n_cells, grid_to_compressed, compressed_to_grid, ...
    nx, ny, dx, dy, cell_area, Z, n_manning, h_0, F_min, ...
    Ks, eta, psi_f, use_infiltration, outlet_direction, slope_outlet, is_outlet_cell)

n = n_cells*2;
res = zeros(n,1);
heps = 1e-9;

for c = 1:n_cells
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1)/nx) + 1;
        j = mod(c-1, nx) + 1;
    end

    b2 = (c-1)*2;
    h  = xHF_new(b2+1);
    F  = xHF_new(b2+2);

    h_old = xHF(b2+1);
    F_old = xHF(b2+2);

    if bdf_order==2 && step>1
        h_prev = xHF_prev(b2+1);
        F_prev = xHF_prev(b2+2);
    else
        h_prev = 0; F_prev = 0;
    end

    [hL, idxL] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'left',  c);
    [hR, idxR] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'right', c);
    [hU, idxU] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'up',    c);
    [hD, idxD] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, i, j, 'down',  c);

    [f_infil, ~, ~] = infiltration_smooth(h, F, dt, i, j, Ks, eta, psi_f, F_min, use_infiltration);

    is_outlet = is_outlet_cell(c);
    nloc = get_manning(n_manning, i, j);

    [Q, ~, ~] = manning_Q_and_dQ(h, hL, hR, hU, hD, Z, i, j, nx, ny, dx, dy, nloc, h_0, ...
        outlet_direction, slope_outlet, is_outlet, heps);

    % inflow (same approximate neighbor calls)
    Q_in = 0;
    if idxL>0, Q_in = Q_in + manning_Q_neighbor_into_us(xHF_new, idxL, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromLeft'); end
    if idxR>0, Q_in = Q_in + manning_Q_neighbor_into_us(xHF_new, idxR, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromRight'); end
    if idxU>0, Q_in = Q_in + manning_Q_neighbor_into_us(xHF_new, idxU, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromUp'); end
    if idxD>0, Q_in = Q_in + manning_Q_neighbor_into_us(xHF_new, idxD, c, Z, nx, ny, dx, dy, n_manning, h_0, ...
            outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, 'fromDown'); end

    Q_out = sum(Q);
    Q_net = Q_in - Q_out;

    dh_dt = R(c) - f_infil + Q_net / cell_area;

    if bdf_order == 1
        res_h = h - h_old - dt*beta*dh_dt;
        res_F = F - F_old - dt*beta*f_infil;
    else
        res_h = h - alpha(1)*h_old - alpha(2)*h_prev - dt*beta*dh_dt;
        res_F = F - alpha(1)*F_old - alpha(2)*F_prev - dt*beta*f_infil;
    end

    res(b2+1) = res_h;
    res(b2+2) = res_F;
end

end

%% ========================================================================
%  
%  ========================================================================
function [t_out, x_full_row] = pack_full_state(t, xHF, params, ...
    use_compressed, n_cells, grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, Z, n_manning, h_0, ...
    outlet_direction, slope_outlet, is_outlet_cell)

heps = 1e-9;
x_full = zeros(n_cells*6,1);

for c=1:n_cells
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1)/nx)+1;
        j = mod(c-1,nx)+1;
    end

    b2 = (c-1)*2;
    h  = xHF(b2+1);
    F  = xHF(b2+2);

    [hL, ~] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'left',  c);
    [hR, ~] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'right', c);
    [hU, ~] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'up',    c);
    [hD, ~] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'down',  c);

    nloc = get_manning(n_manning, i, j);
    [Q, ~, ~] = manning_Q_and_dQ(h, hL, hR, hU, hD, Z, i, j, nx, ny, dx, dy, nloc, h_0, ...
        outlet_direction, slope_outlet, is_outlet_cell(c), heps);

    b6=(c-1)*6;
    x_full(b6+1)=max(h,0);
    x_full(b6+2)=Q(1);
    x_full(b6+3)=Q(2);
    x_full(b6+4)=Q(3);
    x_full(b6+5)=Q(4);
    x_full(b6+6)=F;
end

t_out = t;
x_full_row = x_full(:).';
end

%% ========================================================================
%  Helpers
%  ========================================================================

function nloc = get_manning(n_manning, i, j)
if numel(n_manning) > 1
    nloc = n_manning(i,j);
else
    nloc = n_manning;
end
end

function [hnb, idxnb] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, dir, c)
% returns neighbor h and neighbor cell index (compressed index if used)
idxnb = 0;
hnb = 0;

switch dir
    case 'left'
        if j<=1, return; end
        if use_compressed, idxnb = grid_to_compressed(i, j-1); else, idxnb = (i-1)*nx + (j-1); end
    case 'right'
        if j>=nx, return; end
        if use_compressed, idxnb = grid_to_compressed(i, j+1); else, idxnb = (i-1)*nx + (j+1); end
    case 'up'
        if i<=1, return; end
        if use_compressed, idxnb = grid_to_compressed(i-1, j); else, idxnb = (i-2)*nx + j; end
    case 'down'
        if i>=ny, return; end
        if use_compressed, idxnb = grid_to_compressed(i+1, j); else, idxnb = i*nx + j; end
end

if idxnb>0
    hnb = xHF((idxnb-1)*2 + 1);
end
end

function [f_infil, df_dh, df_dF] = infiltration_smooth(h, F, dt, i, j, Ks, eta, psi_f, F_min, use_infiltration)
% Your same smooth-min infiltration logic, but without hard clamps in residual.
if ~use_infiltration || h <= 0
    f_infil = 0; df_dh=0; df_dF=0; return;
end

% local params
if numel(Ks) > 1
    Ks_c = Ks(i,j);
    eta_c = eta(i,j);
    psi_c = psi_f(i,j);
else
    Ks_c = Ks;
    eta_c = eta;
    psi_c = psi_f;
end

F_safe = max(F, F_min);      % (kept: F_min is model constraint)
eta_h = eta_c * (h + psi_c);
f_cap = Ks_c * (1 + eta_h / F_safe);

i_a = h / dt;

% smooth min: 0.5*(a+b - sqrt((a-b)^2+eps))
% eps must be << (a-b)^2 to avoid corrupting the result
scale = max(f_cap, i_a) + 1e-30;
eps_s = 1e-12 * scale^2;
diff = f_cap - i_a;
s = sqrt(diff*diff + eps_s);
f_infil = 0.5*(f_cap + i_a - s);
f_infil = max(f_infil, 0);

df_dfcap = 0.5*(1 - diff/s);
df_dia   = 0.5*(1 + diff/s);

dfcap_dh = Ks_c * eta_c / F_safe;
dfcap_dF = -Ks_c * eta_h / (F_safe^2);
dia_dh   = 1/dt;

df_dh = df_dfcap*dfcap_dh + df_dia*dia_dh;
df_dF = df_dfcap*dfcap_dF;
end

function [Q, dQdh_own, dQdh_nb] = manning_Q_and_dQ(h, hL, hR, hU, hD, Z, i, j, nx, ny, dx, dy, n, h_0, ...
    outlet_direction, slope_outlet, is_outlet, heps)
% Computes outflows Q=[Ql,Qr,Qu,Qd] from cell (i,j) given own h and neighbor h's
% Returns approximate derivatives wrt own h and each neighbor h (in same order).

% smooth h_eff to avoid singularities in powers/derivs
h_eff = smooth_pos(h - h_0, 1e-6);
h_safe = smooth_pos(h, 1e-3);

z = Z(i,j);
wse = z + h;

% neighbor WSE with boundary/outlet
wseL = wse; wseR = wse; wseU = wse; wseD = wse;

if j>1 && ~isnan(Z(i,j-1)), wseL = Z(i,j-1) + hL;
elseif strcmp(outlet_direction,'left') || is_outlet, wseL = wse - slope_outlet*dx;
end

if j<nx && ~isnan(Z(i,j+1)), wseR = Z(i,j+1) + hR;
elseif strcmp(outlet_direction,'right') || is_outlet, wseR = wse - slope_outlet*dx;
end

if i>1 && ~isnan(Z(i-1,j)), wseU = Z(i-1,j) + hU;
elseif strcmp(outlet_direction,'up') || is_outlet, wseU = wse - slope_outlet*dy;
end

if i<ny && ~isnan(Z(i+1,j)), wseD = Z(i+1,j) + hD;
elseif strcmp(outlet_direction,'down') || is_outlet, wseD = wse - slope_outlet*dy;
end

% Smooth slope positive part
S_l = smooth_pos((wse - wseL)/dx, 1e-8);
S_r = smooth_pos((wse - wseR)/dx, 1e-8);
S_u = smooth_pos((wse - wseU)/dy, 1e-8);
S_d = smooth_pos((wse - wseD)/dy, 1e-8);

S_total = S_l + S_r + S_u + S_d + 1e-12;
S_max = max([S_l,S_r,S_u,S_d]);

v_limit = 8.0;
v_raw = (1/n) * (h_eff^(2/3)) * sqrt(max(S_max,0));
v = min(v_raw, v_limit);

% smooth factor ff (kept from your idea but smoothed)
h_thr = 0.001;
ff = min(1.0, h_safe/h_thr);
ff = ff*ff*(3 - 2*ff);

Q = zeros(4,1);
Q(1) = ff * v * h_eff * dy * (S_l / S_total);
Q(2) = ff * v * h_eff * dy * (S_r / S_total);
Q(3) = ff * v * h_eff * dx * (S_u / S_total);
Q(4) = ff * v * h_eff * dx * (S_d / S_total);

% Approx derivatives:
% own h affects wse and h_eff and v; neighbor h affects slopes (wse_neighbor)
dQdh_own = (5/3) * Q / (h_eff + heps);   % same heuristic but now h_eff is smooth
dQdh_nb  = zeros(4,1);

% neighbor derivative heuristic through slope: when neighbor h increases, slope decreases => Q decreases
% (use your conservative -Q/(S*dx) idea with smoothing)
if S_l > 1e-10, dQdh_nb(1) = -Q(1) / (S_l*dx); end
if S_r > 1e-10, dQdh_nb(2) = -Q(2) / (S_r*dx); end
if S_u > 1e-10, dQdh_nb(3) = -Q(3) / (S_u*dy); end
if S_d > 1e-10, dQdh_nb(4) = -Q(4) / (S_d*dy); end
end

function y = smooth_pos(x, epsv)
% smooth approximation of max(x,0): 0.5*(x + sqrt(x^2 + eps))
y = 0.5*(x + sqrt(x.*x + epsv));
end

function varargout = manning_Q_neighbor_into_us(xHF_new, idxN, idxC, Z, nx, ny, dx, dy, n_manning, h_0, ...
    outlet_direction, slope_outlet, is_outlet_cell, use_compressed, compressed_to_grid, grid_to_compressed, heps, mode)
% Returns neighbor's flux INTO current cell (scalar).
% Also returns approximate dQ/dh_neighbor_own and dQ/dh_current (for Jacobian)
% Usage patterns:
%   Qn = manning_Q_neighbor_into_us(...)                % residual-only calls
%   [Qn, dOwn, dNb] = manning_Q_neighbor_into_us(...)  % Jacobian calls

% Map neighbor idx -> (i,j)
if use_compressed
    in = compressed_to_grid(idxN,1);
    jn = compressed_to_grid(idxN,2);
else
    in = floor((idxN-1)/nx)+1;
    jn = mod((idxN-1),nx)+1;
end

hN = xHF_new((idxN-1)*2 + 1);

% Current cell (receiver) indices for using as "neighbor h" in neighbor's Manning
if use_compressed
    ic = compressed_to_grid(idxC,1);
    jc = compressed_to_grid(idxC,2);
else
    ic = floor((idxC-1)/nx)+1;
    jc = mod((idxC-1),nx)+1;
end
hC = xHF_new((idxC-1)*2 + 1);

% Determine neighbor's four neighbor h's for its Manning call
[hL, ~] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, in, jn, 'left',  idxN);
[hR, ~] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, in, jn, 'right', idxN);
[hU, ~] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, in, jn, 'up',    idxN);
[hD, ~] = neighbor_h(xHF_new, use_compressed, grid_to_compressed, nx, ny, in, jn, 'down',  idxN);

% Override the appropriate neighbor h to be the receiver cell hC (so flux points to us)
switch mode
    case 'fromLeft'   % left neighbor flows right into us -> neighbor's "right neighbor" is us
        hR = hC;
    case 'fromRight'  % right neighbor flows left into us
        hL = hC;
    case 'fromUp'     % up neighbor flows down into us
        hD = hC;
    case 'fromDown'   % down neighbor flows up into us
        hU = hC;
end

% Manning n at neighbor
if numel(n_manning)>1, nloc = n_manning(in,jn); else, nloc = n_manning; end

[Qall, dOwnAll, dNbAll] = manning_Q_and_dQ(hN, hL, hR, hU, hD, Z, in, jn, nx, ny, dx, dy, nloc, h_0, ...
    outlet_direction, slope_outlet, is_outlet_cell(idxN), heps);

% pick component that goes into current cell
switch mode
    case 'fromLeft',  q = Qall(2); dOwn = dOwnAll(2); dNb = dNbAll(2); % neighbor Qr depends on its own h (dOwn) and our h (dNb)
    case 'fromRight', q = Qall(1); dOwn = dOwnAll(1); dNb = dNbAll(1);
    case 'fromUp',    q = Qall(4); dOwn = dOwnAll(4); dNb = dNbAll(4);
    case 'fromDown',  q = Qall(3); dOwn = dOwnAll(3); dNb = dNbAll(3);
end

if nargout <= 1
    varargout{1} = q;
else
    varargout{1} = q;
    varargout{2} = dOwn;
    varargout{3} = dNb;
end

end

function Q_boundary = compute_boundary_outflow_newton(xHF, params, use_compressed, n_cells, ...
    grid_to_compressed, compressed_to_grid, nx, ny, dx, dy, Z, n_manning, h_0, ...
    outlet_direction, slope_outlet, is_outlet_cell)
% Compute total boundary outflow from reduced state [h, F]
% Must reconstruct Q from h using Manning equation
% CRITICAL: Must match boundary treatment in manning_Q_and_dQ exactly

Q_boundary = 0;
heps = 1e-9;

for c = 1:n_cells
    % Get grid position
    if use_compressed
        i = compressed_to_grid(c, 1);
        j = compressed_to_grid(c, 2);
    else
        i = floor((c-1)/nx) + 1;
        j = mod(c-1, nx) + 1;
    end
    
    b2 = (c-1)*2;
    h = xHF(b2+1);
    
    % Get neighbor depths AND check if neighbor exists
    [hL, idxL] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'left',  c);
    [hR, idxR] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'right', c);
    [hU, idxU] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'up',    c);
    [hD, idxD] = neighbor_h(xHF, use_compressed, grid_to_compressed, nx, ny, i, j, 'down',  c);
    
    % Compute Manning Q for this cell
    % Note: manning_Q_and_dQ checks boundary conditions internally
    nloc = get_manning(n_manning, i, j);
    [Q, ~, ~] = manning_Q_and_dQ_boundary_aware(h, hL, hR, hU, hD, idxL, idxR, idxU, idxD, ...
        Z, i, j, nx, ny, dx, dy, nloc, h_0, ...
        outlet_direction, slope_outlet, is_outlet_cell(c), heps);
    
    % Add flows leaving domain boundaries
    % Only count flow where NO valid neighbor exists
    if idxL == 0, Q_boundary = Q_boundary + Q(1); end  % No left neighbor
    if idxR == 0, Q_boundary = Q_boundary + Q(2); end  % No right neighbor
    if idxU == 0, Q_boundary = Q_boundary + Q(3); end  % No up neighbor
    if idxD == 0, Q_boundary = Q_boundary + Q(4); end  % No down neighbor
end

end

function [Q, dQdh_own, dQdh_nb] = manning_Q_and_dQ_boundary_aware(h, hL, hR, hU, hD, idxL, idxR, idxU, idxD, ...
    Z, i, j, nx, ny, dx, dy, n, h_0, outlet_direction, slope_outlet, is_outlet, heps)
% Computes outflows Q=[Ql,Qr,Qu,Qd] with proper boundary condition handling
% Takes neighbor indices to determine if neighbor exists (idx>0) or is boundary (idx==0)

h_eff = smooth_pos(h - h_0, 1e-6);
h_safe = smooth_pos(h, 1e-3);

z = Z(i,j);
wse = z + h;

% Neighbor WSE: use actual neighbor if exists, else apply boundary condition
wseL = wse; wseR = wse; wseU = wse; wseD = wse;

% Left
if idxL > 0
    % Valid neighbor exists
    if ~isnan(Z(i,j-1))
        wseL = Z(i,j-1) + hL;
    end
elseif strcmp(outlet_direction,'left') || is_outlet
    % Boundary with outlet
    wseL = wse - slope_outlet*dx;
elseif j > 1 && ~isnan(Z(i,j-1))
    % Boundary but terrain exists - use terrain slope for free drainage
    terrain_slope = max((z - Z(i,j-1))/dx, 0.001);
    wseL = wse - terrain_slope*dx;
end

% Right
if idxR > 0
    if ~isnan(Z(i,j+1))
        wseR = Z(i,j+1) + hR;
    end
elseif strcmp(outlet_direction,'right') || is_outlet
    wseR = wse - slope_outlet*dx;
elseif j < nx && ~isnan(Z(i,j+1))
    terrain_slope = max((z - Z(i,j+1))/dx, 0.001);
    wseR = wse - terrain_slope*dx;
end

% Up
if idxU > 0
    if ~isnan(Z(i-1,j))
        wseU = Z(i-1,j) + hU;
    end
elseif strcmp(outlet_direction,'up') || is_outlet
    wseU = wse - slope_outlet*dy;
elseif i > 1 && ~isnan(Z(i-1,j))
    terrain_slope = max((z - Z(i-1,j))/dy, 0.001);
    wseU = wse - terrain_slope*dy;
end

% Down
if idxD > 0
    if ~isnan(Z(i+1,j))
        wseD = Z(i+1,j) + hD;
    end
elseif strcmp(outlet_direction,'down') || is_outlet
    wseD = wse - slope_outlet*dy;
elseif i < ny && ~isnan(Z(i+1,j))
    terrain_slope = max((z - Z(i+1,j))/dy, 0.001);
    wseD = wse - terrain_slope*dy;
end

% Rest is same as manning_Q_and_dQ
S_l = smooth_pos((wse - wseL)/dx, 1e-8);
S_r = smooth_pos((wse - wseR)/dx, 1e-8);
S_u = smooth_pos((wse - wseU)/dy, 1e-8);
S_d = smooth_pos((wse - wseD)/dy, 1e-8);

S_total = S_l + S_r + S_u + S_d + 1e-12;
S_max = max([S_l,S_r,S_u,S_d]);

v_limit = 8.0;
v_raw = (1/n) * (h_eff^(2/3)) * sqrt(max(S_max,0));
v = min(v_raw, v_limit);

h_thr = 0.001;
ff = min(1.0, h_safe/h_thr);
ff = ff*ff*(3 - 2*ff);

Q = zeros(4,1);
Q(1) = ff * v * h_eff * dy * (S_l / S_total);
Q(2) = ff * v * h_eff * dy * (S_r / S_total);
Q(3) = ff * v * h_eff * dx * (S_u / S_total);
Q(4) = ff * v * h_eff * dx * (S_d / S_total);

dQdh_own = (5/3) * Q / (h_eff + heps);
dQdh_nb  = zeros(4,1);

if S_l > 1e-10, dQdh_nb(1) = -Q(1) / (S_l*dx); end
if S_r > 1e-10, dQdh_nb(2) = -Q(2) / (S_r*dx); end
if S_u > 1e-10, dQdh_nb(3) = -Q(3) / (S_u*dy); end
if S_d > 1e-10, dQdh_nb(4) = -Q(4) / (S_d*dy); end
end

