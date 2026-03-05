function params = add_infiltration_params(params, Ks, eta, psi_f, F_min)
% add_infiltration_params: Add Green-Ampt infiltration parameters to params
%
% This function adds infiltration parameters to an existing params structure
% for use with the 6-state DAE formulation (hydrodynamic + hydrologic).
%
% Inputs:
%   params - Existing parameters structure (from domain setup)
%   Ks - Saturated hydraulic conductivity [m/s]
%        Can also be specified as [mm/h] and will be converted
%   eta - Effective porosity minus initial moisture content [-]
%         Typical range: 0.2-0.5
%   psi_f - Wetting front suction head [m]
%           Typical range: 0.05-0.3 m
%   F_min - Minimum cumulative infiltration to prevent stiffness [m]
%           Recommended: 0.010 m (10 mm)
%
% Author: M. H. Kazma
% Date: November 8, 2025

%% Input validation
if nargin < 5
    error('All 5 arguments required: params, Ks, eta, psi_f, F_min');
end

%% Auto-detect if Ks is in mm/h and convert to m/s
if Ks > 0.01  % Likely in mm/h if > 0.01
    fprintf('Note: Ks appears to be in mm/h, converting to m/s\n');
    Ks_mms = Ks / 3600 / 1000;  % mm/h -> m/s
    fprintf('  Ks = %.2f mm/h = %.2e m/s\n', Ks, Ks_mms);
    Ks = Ks_mms;
end

%% Add infiltration parameters
params.Ks = Ks;
params.eta = eta;
params.psi_f = psi_f;
params.F_min = F_min;
if ~isfield(params, 'F_initial')
    params.F_initial = F_min;  % Default: start at F_min (dry antecedent)
end

%% Display summary
fprintf('\nInfiltration parameters added:\n');
fprintf('  Ks = %.2e m/s (%.2f mm/h)\n', Ks, Ks*3600*1000);
fprintf('  eta = %.3f (effective porosity)\n', eta);
fprintf('  psi_f = %.3f m (wetting front suction)\n', psi_f);
fprintf('  F_min = %.3f m (%.1f mm) - prevents stiffness\n\n', F_min, F_min*1000);

% Provide guidance on soil type
Ks_mmh = Ks * 3600 * 1000;
if Ks_mmh > 10
    soil_type = 'Sand/Sandy loam (high permeability)';
elseif Ks_mmh > 1
    soil_type = 'Loam/Silt loam (moderate permeability)';
else
    soil_type = 'Clay/Clay loam (low permeability)';
end
fprintf('Approximate soil type: %s\n\n', soil_type);

end
