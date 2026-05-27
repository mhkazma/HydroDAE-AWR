"""
Post-processing for AWR2026 R1 Reviewer Response, Comment C1-R2 (cell-track).

Companion to AWR2026_postprocess_C1_R2.py.  Reads the persisted N_MC = 100
tracked-cell MC ensemble of depth h(t) at representative gauged and
ungauged interior cells (the same cells visualized in Figs. 11-12), and
computes the phase-stratified coverage of the proposed Gaussian-rendered
covariance interval mu_h(t) +/- 1.96 sigma_h(t) cell-by-cell.

Phases match AWR2026_postprocess_C1_R2.py:
    rising     : t <  t_peak     and  mu_h(t) <  (1 - w) * h_peak
    peak       : mu_h(t)       >= (1 - w) * h_peak     (w = 0.10)
    recession  : t >  t_peak     and  mu_h(t) <  (1 - w) * h_peak

Inputs (already on disk; no re-simulation):
    results/<case>/AWR2026/<case>_mc_celltrack_ca_routing_ALL.mat

Output (per watershed):
    results/<case>/AWR2026/<case>_R1_C1R2_celltrack_coverage.mat
"""

import os
import numpy as np
import h5py
from scipy.io import savemat


def _read_str_attr(h, key):
    """Read a MATLAB string/char array from an HDF5 dataset path."""
    try:
        d = h[key]
    except KeyError:
        return None
    arr = np.array(d)
    if arr.dtype.kind in ('U', 'S'):
        return ''.join(arr.flatten().tolist())
    # uint16 char array
    try:
        return ''.join(chr(int(c)) for c in arr.flatten() if int(c) != 0)
    except Exception:
        return None


def _read_track_labels(h):
    """Read mc_results.track_labels (MATLAB cellstr) via HDF5 references."""
    try:
        refs = h['mc_results/track_labels']
    except KeyError:
        return None
    labels = []
    for r in np.array(refs).flatten():
        try:
            obj = h[r]
            txt = ''.join(chr(int(c)) for c in np.array(obj).flatten() if int(c) != 0)
            labels.append(txt)
        except Exception:
            labels.append('')
    return labels


def _load_celltrack(path):
    """
    Load the cell-track MC ensemble.

    MATLAB save shapes:
      mc_results.h_track_ensemble : [N_ok, n_track, nt]  -> HDF5 [nt, n_track, N_ok]
      mc_results.h_track_mean     : [n_track, nt]        -> HDF5 [nt, n_track]
      mc_results.h_track_std      : [n_track, nt]        -> HDF5 [nt, n_track]
      mc_results.t_common         : [nt, 1]              -> HDF5 [1, nt]
      mc_results.track_cells      : [n_track, 2]         -> HDF5 [2, n_track]
    """
    with h5py.File(path, 'r') as h:
        H_ens = np.array(h['mc_results/h_track_ensemble'])  # [nt, n_track, N]
        H_mean = np.array(h['mc_results/h_track_mean'])      # [nt, n_track]
        H_std = np.array(h['mc_results/h_track_std'])        # [nt, n_track]
        t = np.array(h['mc_results/t_common']).squeeze()
        tc = np.array(h['mc_results/track_cells'])           # [2, n_track]
        labels = _read_track_labels(h)

    # Re-orient to Pythonic [N, n_track, nt] / [n_track, nt] / [n_track, 2]
    if H_ens.ndim != 3:
        raise RuntimeError(f'h_track_ensemble has ndim={H_ens.ndim}')
    nt = t.size
    # HDF5 stores MATLAB [N, n_track, nt] in reversed order [nt, n_track, N].
    if H_ens.shape[0] == nt:
        H_ens = np.transpose(H_ens, (2, 1, 0))  # -> [N, n_track, nt]
    if H_mean.shape[0] == nt:
        H_mean = H_mean.T                       # -> [n_track, nt]
        H_std = H_std.T
    if tc.shape[0] == 2:
        tc = tc.T                               # -> [n_track, 2]

    return t, H_ens, H_mean, H_std, tc.astype(int), labels


def cell_coverage_by_phase(t, h_ens, mu, sd, peak_window_frac=0.10, k=1.96):
    """
    Phase-stratified coverage for one tracked cell.

    h_ens : [N, nt]   ensemble depth at the cell
    mu    : [nt]      MC mean
    sd    : [nt]      MC std
    """
    N, nt = h_ens.shape
    ipk = int(np.argmax(mu))
    h_pk = float(mu[ipk])

    lo = mu - k * sd
    hi = mu + k * sd
    inside = (h_ens >= lo[None, :]) & (h_ens <= hi[None, :])  # [N, nt]

    if h_pk > 0:
        near_peak = mu >= (1.0 - peak_window_frac) * h_pk
    else:
        near_peak = np.zeros_like(mu, dtype=bool)
    rising = (t < t[ipk]) & ~near_peak
    peak = near_peak
    recess = (t > t[ipk]) & ~near_peak

    def _stat(mask):
        if not np.any(mask):
            return dict(n_t=0, n_pairs=0, coverage=float('nan'),
                        t_start=float('nan'), t_end=float('nan'))
        block = inside[:, mask]
        return dict(n_t=int(mask.sum()),
                    n_pairs=int(block.size),
                    coverage=float(block.mean()),
                    t_start=float(t[mask].min()),
                    t_end=float(t[mask].max()))

    return dict(
        n_MC=int(N), T=int(nt),
        t_peak=float(t[ipk]), h_peak=h_pk, peak_index=int(ipk),
        peak_window_frac=float(peak_window_frac), k_sigma=float(k),
        rising=_stat(rising),
        peak=_stat(peak),
        recession=_stat(recess),
        overall=_stat(np.ones_like(mu, dtype=bool)),
    )


ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
CASES = [
    ('V-Tilted',     'vtilted',      'vtilted'),
    ('Walnut Gulch', 'walnut_gulch', 'walnut_gulch'),
]


def _paths(case_dir, file_prefix):
    base = os.path.join(ROOT, 'results', case_dir, 'AWR2026')
    return dict(
        mc=os.path.join(base, f'{file_prefix}_mc_celltrack_ca_routing_ALL.mat'),
        out=os.path.join(base, f'{file_prefix}_R1_C1R2_celltrack_coverage.mat'),
    )


def main():
    for name, case_dir, prefix in CASES:
        paths = _paths(case_dir, prefix)
        print('=' * 72)
        print(f' [{name}]')
        print('=' * 72)
        if not os.path.isfile(paths['mc']):
            print(f'  MISSING: {paths["mc"]} -- skipping')
            continue
        t, H_ens, H_mean, H_std, tc, labels = _load_celltrack(paths['mc'])
        n_track = tc.shape[0]
        if labels is None or len(labels) != n_track:
            labels = [f'cell_{i+1}' for i in range(n_track)]

        per_cell = []
        for j in range(n_track):
            cov = cell_coverage_by_phase(t, H_ens[:, j, :], H_mean[j, :], H_std[j, :])
            cov['label'] = labels[j]
            cov['row'] = int(tc[j, 0])
            cov['col'] = int(tc[j, 1])
            per_cell.append(cov)

        # Console print
        hdr = f'  {"label":<12s} (row,col)   rising   peak   recession   overall'
        print(hdr)
        print('  ' + '-' * (len(hdr) - 2))
        for c in per_cell:
            def _pct(x):
                return f'{100.0*x:5.1f}%' if not np.isnan(x) else '  N/A '
            print(f'  {c["label"]:<12s} ({c["row"]:3d},{c["col"]:3d}) '
                  f'  {_pct(c["rising"]["coverage"])} '
                  f' {_pct(c["peak"]["coverage"])} '
                  f'   {_pct(c["recession"]["coverage"])} '
                  f'   {_pct(c["overall"]["coverage"])}')

        mdict = {
            'per_cell_coverage': per_cell,
            'track_cells': tc,
            'track_labels': np.array(labels, dtype=object),
            'meta': dict(
                response='AWR2026 R1 -- C1-R2 cell-track coverage',
                case=name,
                sources=dict(mc_file=os.path.basename(paths['mc'])),
                peak_window_frac=0.10,
                k_sigma=1.96,
                interval='mu_h(t) +/- 1.96 sigma_h(t) per cell (Gaussian-'
                         'rendered covariance interval from the LHS MC '
                         'ensemble at representative gauged and ungauged '
                         'interior cells from Figs. 11-12).',
            ),
        }
        savemat(paths['out'], mdict, do_compression=True, oned_as='column')
        print(f'  wrote {paths["out"]}')


if __name__ == '__main__':
    main()
