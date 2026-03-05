function Z_filled = fill_dem_pits(Z, valid_mask)
%FILL_DEM_PITS Fill interior pits/depressions in a DEM.
%
%   Z_filled = fill_dem_pits(Z, valid_mask)
%
%   Raises the elevation of interior pit cells to the minimum spill level
%   so that every valid cell has a non-ascending path to the watershed
%   boundary. Uses iterative flooding.
%
%   Inputs:
%       Z          - [nr x nc] DEM (NaN outside watershed)
%       valid_mask - [nr x nc] logical (true = inside watershed)
%
%   Output:
%       Z_filled   - [nr x nc] pit-filled DEM

[nr, nc] = size(Z);
eps_fill = 1e-4;  % 0.1 mm increment to ensure strict drainage

% Neighbor offsets (4-connected)
di = [-1, 1, 0, 0];
dj = [0, 0, -1, 1];

% ---- Step 1: Initialize water surface ----
% Boundary cells of the valid domain: set to their own elevation.
% All other valid cells: set to infinity (will be lowered during flooding).
INF_ELEV = 1e20;
W = INF_ELEV * ones(nr, nc);

is_boundary = false(nr, nc);
for i = 1:nr
    for j = 1:nc
        if ~valid_mask(i,j), continue; end
        
        on_edge = (i == 1 || i == nr || j == 1 || j == nc);
        if ~on_edge
            for d = 1:4
                ni = i + di(d);
                nj = j + dj(d);
                if ~valid_mask(ni, nj)
                    on_edge = true;
                    break;
                end
            end
        end
        
        if on_edge
            W(i,j) = Z(i,j);
            is_boundary(i,j) = true;
        end
    end
end

n_boundary = sum(is_boundary(:));

% ---- Step 2: Iterative flooding (forward + backward sweeps) ----
changed = true;
n_iter = 0;

while changed
    changed = false;
    n_iter = n_iter + 1;
    
    % Forward sweep
    for i = 1:nr
        for j = 1:nc
            if ~valid_mask(i,j) || is_boundary(i,j), continue; end
            
            for d = 1:4
                ni = i + di(d);
                nj = j + dj(d);
                if ni < 1 || ni > nr || nj < 1 || nj > nc, continue; end
                if ~valid_mask(ni, nj), continue; end
                
                if Z(i,j) >= W(ni,nj) + eps_fill
                    % Cell is naturally above neighbor's water level
                    if W(i,j) > Z(i,j)
                        W(i,j) = Z(i,j);
                        changed = true;
                    end
                else
                    % Cell is in a depression: set to neighbor water + eps
                    new_w = W(ni,nj) + eps_fill;
                    if new_w < W(i,j)
                        W(i,j) = new_w;
                        changed = true;
                    end
                end
            end
        end
    end
    
    % Backward sweep for faster convergence
    for i = nr:-1:1
        for j = nc:-1:1
            if ~valid_mask(i,j) || is_boundary(i,j), continue; end
            
            for d = 1:4
                ni = i + di(d);
                nj = j + dj(d);
                if ni < 1 || ni > nr || nj < 1 || nj > nc, continue; end
                if ~valid_mask(ni, nj), continue; end
                
                if Z(i,j) >= W(ni,nj) + eps_fill
                    if W(i,j) > Z(i,j)
                        W(i,j) = Z(i,j);
                        changed = true;
                    end
                else
                    new_w = W(ni,nj) + eps_fill;
                    if new_w < W(i,j)
                        W(i,j) = new_w;
                        changed = true;
                    end
                end
            end
        end
    end
end

% ---- Step 3: Apply fill and compute stats ----
Z_filled = Z;
n_raised = 0;
max_raise = 0;
sum_raise = 0;

for i = 1:nr
    for j = 1:nc
        if valid_mask(i,j) && W(i,j) < INF_ELEV
            raise = W(i,j) - Z(i,j);
            if raise > eps_fill
                n_raised = n_raised + 1;
                max_raise = max(max_raise, raise);
                sum_raise = sum_raise + raise;
            end
            Z_filled(i,j) = W(i,j);
        end
    end
end

n_valid = sum(valid_mask(:));
if n_raised > 0
    mean_raise = sum_raise / n_raised;
else
    mean_raise = 0;
end

fprintf('  Pit-fill: %d boundary seeds, %d iterations\n', n_boundary, n_iter);
fprintf('  %d cells raised (%.1f%%), max raise=%.2f m, mean raise=%.2f m\n', ...
    n_raised, 100*n_raised/n_valid, max_raise, mean_raise);

end
