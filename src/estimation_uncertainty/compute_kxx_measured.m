function [K_xx, K_yy, rc, had_neg_diag] = compute_kxx_measured(A_dae, K_bb, C_local, K_vv_local, W)
% COMPUTE_KXX_MEASURED  PSE covariance with measurement augmentation
%
% USAGE:
%   [K_xx, K_yy, rc, had_neg_diag] = compute_kxx_measured(A_dae, K_bb, C_local, K_vv_local, W)
%
% DESCRIPTION:
%   Implements the PSE Theorem 1 with measurement
%   augmentation. The system matrix A^s is formed by stacking the DAE
%   Jacobian E with the measurement selection matrix C:
%
%       A^s = [E; C]    (over-determined if C is non-empty)
%
%   The covariance of the RHS becomes:
%
%       K_bb^s = blkdiag(K_bb, K_vv)
%
%   And the GLS/BLUE solution (Corollary 1 with W = (K_bb^s)^{-1}):
%
%       K_xx = (A^T W A)^{-1}    (guaranteed variance reduction)
%
%   For the over-determined case with custom weighting (Corollary 1):
%
%       K_xx = (A^T W A)^{-1} A^T W K_bb^s W A (A^T W A)^{-1}
%
%   The output measurement covariance is:
%
%       K_yy = C * K_xx * C' + K_vv
%
% MATHEMATICAL INSIGHT:
%   Adding measurements (rows in C) to the system matrix A^s makes it
%   over-determined (more equations than unknowns). This REDUCES the
%   state covariance K_xx because additional information constrains the
%   state estimates. The noise variance K_vv represents measurement
%   imprecision — better sensors (smaller K_vv) yield tighter bounds.
%
% INPUTS:
%   A_dae      - n_x × n_x DAE Jacobian (system matrix E)
%   K_bb       - n_x × n_x covariance of DAE RHS (parameter + temporal)
%   C_local    - n_y × n_x measurement matrix ([] if no sensors)
%   K_vv_local - n_y × n_y measurement noise covariance ([] if no sensors)
%   W          - Weight matrix for over-determined case ([] for unweighted)
%
% OUTPUTS:
%   K_xx       - n_x × n_x state covariance matrix
%   K_yy       - n_y × n_y output (measurement) covariance matrix
%                ([] if no sensors)
%   rc         - Reciprocal condition number of system matrix
%   had_neg_diag - True if eigenvalue cleanup was needed
%
% REFERENCE:
%   Theorem 1 (Eq. 14) and Corollary 1 (Eq. 18)
%
% Author: M. H. Kazma
% Date: February 2026

n_x = size(A_dae, 1);
has_measurements = ~isempty(C_local) && size(C_local, 1) > 0;

if has_measurements
    n_y = size(C_local, 1);
    
    %% Build augmented system: A^s = [E; C], K_bb^s = blkdiag(K_bb, K_vv)
    % Equation (10):
    %   A^s(k) = [E(k); C]
    %   b^s(k) = [z(k); y(k) - v(k)]
    A_aug = [A_dae; C_local];          % (n_x + n_y) × n_x
    K_bb_aug = blkdiag(K_bb, K_vv_local);  % (n_x + n_y) × (n_x + n_y)
    
    % Ensure symmetry
    K_bb_aug = 0.5 * (K_bb_aug + K_bb_aug');
    
    %% Compute GLS weight matrix W = (K_bb^s)^{-1}
    % The BLUE/Aitken estimator for the over-determined system:
    %   K_xx = (A' W A)^{-1}   where W = (K_bb^s)^{-1}
    % This GUARANTEES variance reduction: adding measurements can only
    % decrease K_xx in the positive-semidefinite sense.
    %
    % K_bb^s = blkdiag(K_bb, K_vv), so:
    %   W = blkdiag(inv(K_bb), inv(K_vv))
    
    if isempty(W)
        % Compute proper block-diagonal inverse
        % K_bb is 6x6 (cheap to invert), K_vv is diagonal
        rc_bb = rcond(K_bb);
        if rc_bb > 1e-14 && ~isnan(rc_bb)
            K_bb_inv = K_bb \ eye(n_x);
        else
            % Regularize singular K_bb
            lambda = 1e-10 * max(max(abs(diag(K_bb))), 1e-20);
            K_bb_inv = (K_bb + lambda * eye(n_x)) \ eye(n_x);
        end
        K_vv_inv = diag(1 ./ max(diag(K_vv_local), 1e-30));
        W_full = blkdiag(K_bb_inv, K_vv_inv);
    else
        % User-provided weight matrix
        n_aug = size(A_aug, 1);
        if isscalar(W)
            W_full = W * eye(n_aug);
        elseif numel(W) == n_aug
            W_full = diag(W(:));
        elseif all(size(W) == n_aug)
            W_full = W;
        else
            % Fallback: auto-GLS
            rc_bb = rcond(K_bb);
            if rc_bb > 1e-14, K_bb_inv = K_bb \ eye(n_x);
            else, K_bb_inv = (K_bb + 1e-10*eye(n_x)) \ eye(n_x); end
            K_vv_inv = diag(1 ./ max(diag(K_vv_local), 1e-30));
            W_full = blkdiag(K_bb_inv, K_vv_inv);
        end
    end
    
    %% GLS solution: K_xx = (A' W A)^{-1}  (BLUE estimator)
    % When W = (K_bb^s)^{-1}, the sandwich formula collapses:
    %   (A'WA)^{-1} A'W K_bb W A (A'WA)^{-1} = (A'WA)^{-1}
    % So we compute K_xx = inv(A' W A) directly.
    A_W = A_aug' * W_full * A_aug;   % n_x × n_x  (Fisher information)
    rc = rcond(A_W);
    
    if rc > 1e-12 && ~isnan(rc)
        K_xx = A_W \ eye(n_x);
    else
        % Tikhonov regularization
        lambda = 1e-10 * max(max(abs(diag(A_W))), 1e-16);
        A_W_reg = A_W + lambda * eye(n_x);
        K_xx = A_W_reg \ eye(n_x);
    end
    
    % Ensure symmetry
    K_xx = 0.5 * (K_xx + K_xx');
    
    %% Output measurement covariance: K_yy = C K_xx C' + K_vv
    K_yy = C_local * K_xx * C_local' + K_vv_local;
    K_yy = 0.5 * (K_yy + K_yy');
    
else
    %% No measurements: fall back to standard PSE
    % K_xx = A^{-1} K_bb A^{-T}  (square system)
    K_bb = 0.5 * (K_bb + K_bb');
    rc = rcond(A_dae);
    
    if rc > 1e-12 && ~isnan(rc)
        A_inv = A_dae \ eye(n_x);
        K_xx = A_inv * K_bb * A_inv';
    else
        AtA = A_dae' * A_dae;
        lambda = 1e-10 * max(max(abs(diag(AtA))), 1e-16);
        AtA_reg = AtA + lambda * eye(n_x);
        AtA_inv = AtA_reg \ eye(n_x);
        AtA_inv_At = AtA_inv * A_dae';
        K_xx = AtA_inv_At * K_bb * AtA_inv_At';
    end
    
    K_xx = 0.5 * (K_xx + K_xx');
    K_yy = [];
end

%% Eigenvalue cleanup: clamp negative diagonal (numerical artifact)
had_neg_diag = any(diag(K_xx) < 0);
if had_neg_diag
    [V, D] = eig(K_xx);
    D = diag(max(diag(D), 0));
    K_xx = V * D * V';
    K_xx = 0.5 * (K_xx + K_xx');
end

end
