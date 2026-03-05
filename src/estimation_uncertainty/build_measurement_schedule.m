function schedule = build_measurement_schedule(t_out, varargin)
% BUILD_MEASUREMENT_SCHEDULE  Generate temporal measurement schedule for PSE
%
% USAGE:
%   schedule = build_measurement_schedule(t_out, 'Interval', 60)
%   schedule = build_measurement_schedule(t_out, 'Times', [100, 200, 300])
%   schedule = build_measurement_schedule(t_out, 'Fraction', 0.1)
%
% DESCRIPTION:
%   Creates a logical schedule vector indicating which timesteps have
%   measurement data available. This enables the "sawtooth" variance
%   pattern in PSE: at measurement times, the augmented system
%   A^s = [E; C] is used (variance drops), while at non-measurement
%   times only the physics A^s = E propagates uncertainty (variance grows).
%
%   This temporal scheduling is the key scientific contribution of the
%   AWR2026 paper — it captures the realistic scenario where sensor
%   measurements arrive at discrete intervals, not continuously.
%
% INPUTS:
%   t_out    - Time vector from simulation [nt x 1] or [1 x nt], in seconds
%
% NAME-VALUE PAIRS:
%   'Interval'  - Measurement interval in seconds (e.g., 60 = every minute)
%                 Measurements occur at t = 0, Interval, 2*Interval, ...
%   'Times'     - Explicit measurement times in seconds [1 x n_meas]
%                 Each time is snapped to the nearest simulation timestep
%   'Fraction'  - Fraction of timesteps with measurements (0 < f <= 1)
%                 E.g., 0.1 means ~10% of timesteps have measurements
%                 Timesteps are uniformly spaced
%   'Indices'   - Direct timestep indices where measurements are available
%                 [1 x n_meas] integer indices into t_out
%   'Always'    - If true, measurements at every timestep (default: false)
%                 This reproduces the original (non-scheduled) behavior
%
% OUTPUTS:
%   schedule    - struct with fields:
%     .is_measurement  - [1 x nt] logical vector
%     .measurement_indices - sorted vector of measurement timestep indices
%     .n_measurements  - total number of measurement timesteps
%     .n_total         - total timesteps
%     .fraction        - fraction of timesteps with measurements
%     .description     - human-readable description of the schedule
%
% EXAMPLES:
%   % Measurement every 60 seconds over a 3600-second simulation
%   t = linspace(0, 3600, 7200);
%   s = build_measurement_schedule(t, 'Interval', 60);
%   % s.n_measurements ≈ 61, s.fraction ≈ 0.85%
%
%   % 10% of timesteps have measurements
%   s = build_measurement_schedule(t, 'Fraction', 0.1);
%
% REFERENCE:
%   AWR2026 paper, Section IV-B: Temporal measurement scheduling
%
% Author: M. H. Kazma
% Date: 2025

%% Parse inputs
p = inputParser;
p.addRequired('t_out', @isnumeric);
p.addParameter('Interval', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('Times', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('Fraction', [], @(x) isempty(x) || (isscalar(x) && x > 0 && x <= 1));
p.addParameter('Indices', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('Always', false, @islogical);
p.parse(t_out, varargin{:});
opts = p.Results;

nt = numel(t_out);
t_out = t_out(:)';  % ensure row vector

%% Count how many modes are specified
n_modes = ~isempty(opts.Interval) + ~isempty(opts.Times) + ...
          ~isempty(opts.Fraction) + ~isempty(opts.Indices) + opts.Always;

if n_modes == 0
    % Default: always available (backward-compatible)
    opts.Always = true;
end

if n_modes > 1 && ~opts.Always
    warning('build_measurement_schedule:multipleInputs', ...
        'Multiple schedule modes specified. Using first found: Interval > Times > Fraction > Indices.');
end

%% Build schedule
is_measurement = false(1, nt);

if ~isempty(opts.Interval) && ~opts.Always
    % Interval-based: measurements at t = 0, dt_meas, 2*dt_meas, ...
    dt_meas = opts.Interval;
    t_start = t_out(1);
    meas_times = t_start:dt_meas:(t_out(end) + dt_meas/2);
    
    for mi = 1:numel(meas_times)
        [~, idx] = min(abs(t_out - meas_times(mi)));
        is_measurement(idx) = true;
    end
    desc = sprintf('Interval-based: every %.1f s (%d measurements in %.1f s)', ...
        dt_meas, nnz(is_measurement), t_out(end) - t_out(1));

elseif ~isempty(opts.Times) && ~opts.Always
    % Explicit times: snap to nearest simulation timestep
    meas_times = opts.Times(:)';
    for mi = 1:numel(meas_times)
        [~, idx] = min(abs(t_out - meas_times(mi)));
        is_measurement(idx) = true;
    end
    desc = sprintf('Explicit times: %d measurement times specified', numel(meas_times));

elseif ~isempty(opts.Fraction) && ~opts.Always
    % Fraction-based: uniformly spaced subset
    n_meas = max(1, round(nt * opts.Fraction));
    meas_idx = round(linspace(1, nt, n_meas));
    is_measurement(meas_idx) = true;
    desc = sprintf('Fraction-based: %.1f%% of timesteps (%d/%d)', ...
        100 * opts.Fraction, n_meas, nt);

elseif ~isempty(opts.Indices) && ~opts.Always
    % Direct indices
    idx = opts.Indices(:)';
    idx = idx(idx >= 1 & idx <= nt);
    is_measurement(idx) = true;
    desc = sprintf('Index-based: %d measurement timesteps', numel(idx));

else
    % Always: measurements at every timestep
    is_measurement(:) = true;
    desc = 'Continuous: measurements at every timestep (no scheduling)';
end

% Ensure at least first and last timestep have measurements for stability
is_measurement(1) = true;
is_measurement(nt) = true;

%% Package output
schedule = struct();
schedule.is_measurement = is_measurement;
schedule.measurement_indices = find(is_measurement);
schedule.n_measurements = nnz(is_measurement);
schedule.n_total = nt;
schedule.fraction = schedule.n_measurements / nt;
schedule.description = desc;

fprintf('  Measurement schedule: %s\n', desc);
fprintf('  %d/%d timesteps have measurements (%.1f%%)\n', ...
    schedule.n_measurements, nt, 100 * schedule.fraction);

end
