"""
Post-processing for AWR2026 Reviewer Response (R1).

Computes the numerical values needed for:
  - C2-R1: MC sample skewness gamma_1 and excess kurtosis gamma_2 of the
           outlet discharge Q at the peak time, and the empirical coverage
           of the covariance-based mu +/- 1.96 sigma interval against the
           nominal 95% level, for V-Tilted and Walnut Gulch.
  - C3-R1: NSE and KGE of the DAE outlet hydrograph against the baseline
           solvers (CA routing, semi-implicit local inertia) for V-Tilted
           and Walnut Gulch.

Inputs (already on disk; no re-simulation):
  results/vtilted/AWR2026/vtilted_mc_validation_ca_routing_ALL.mat
  results/vtilted/AWR2026/vtilted_pse_measured_ALL.mat
  results/vtilted/AWR2026/vtilted_solver_comparison.mat
  results/walnut_gulch/AWR2026/walnut_gulch_mc_validation_ca_routing_ALL.mat
  results/walnut_gulch/AWR2026/walnut_gulch_pse_measured_ALL.mat
  results/walnut_gulch/AWR2026/walnut_gulch_solver_comparison.mat
"""

import os
import numpy as np
import h5py
from scipy.io import savemat


def read_str(h, ds):
    """Decode a MATLAB v7.3 char array stored as uint16."""
    arr = np.array(ds).squeeze()
    return ''.join(chr(int(c)) for c in arr.flat)


def deref(h, ref):
    """Dereference an HDF5 object reference."""
    return np.array(h[ref])


def load_mc_outlet(path):
    with h5py.File(path, 'r') as h:
        Q_ens = np.array(h['mc_results/Q_outlet_ensemble'])  # (n_MC, T) or (T, n_MC)
        Q_mean = np.array(h['mc_results/Q_outlet_mean']).squeeze()
        Q_std = np.array(h['mc_results/Q_outlet_std']).squeeze()
        t = np.array(h['mc_results/t_common']).squeeze()
    # h5py reads MATLAB (T, n_MC) as (n_MC, T); confirm orientation by matching mean length.
    if Q_ens.shape[0] == Q_mean.size:
        # Q_ens is (T, n_MC) after transpose to MATLAB layout
        Q_ens = Q_ens.T  # now (n_MC, T)
    return t, Q_mean, Q_std, Q_ens


def load_pse_outlet_std(path, t_target):
    """Return the PSE std of outlet Q on the same time grid as MC."""
    with h5py.File(path, 'r') as h:
        # pse_meas/Q_d_std has shape (sensors, grid_y*grid_x?, T) in MATLAB.
        # We use Q_outlet equivalent from y_std if it exists.
        # Simpler: load var_y at the outlet sensor (last sensor by convention).
        # Fall back to ensemble std if PSE-side outlet not directly available.
        t_pse = np.array(h['pse_meas/t']).squeeze()
        # y_std is (sensors, T) in MATLAB -> (T, sensors) after h5py read
        y_std = np.array(h['pse_meas/y_std'])
    return t_pse, y_std


def skewness_kurtosis_coverage(Q_ens_peak, mu, sigma, k=1.96):
    """Sample skewness, excess kurtosis, and empirical coverage of [mu+/-k*sigma]."""
    x = np.asarray(Q_ens_peak, dtype=float)
    n = x.size
    xm = x - x.mean()
    s2 = (xm ** 2).mean()
    s = np.sqrt(s2)
    if s == 0.0:
        return 0.0, 0.0, 1.0
    m3 = (xm ** 3).mean()
    m4 = (xm ** 4).mean()
    g1 = m3 / s ** 3                          # sample skewness (Fisher)
    g2 = m4 / s ** 4 - 3.0                    # excess kurtosis
    lo, hi = mu - k * sigma, mu + k * sigma
    cov = float(((x >= lo) & (x <= hi)).mean())
    return g1, g2, cov


def nse(sim, obs):
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    denom = ((obs - obs.mean()) ** 2).sum()
    if denom == 0.0:
        return float('nan')
    return 1.0 - ((sim - obs) ** 2).sum() / denom


def kge(sim, obs):
    """Kling-Gupta Efficiency (Gupta et al. 2009)."""
    sim, obs = np.asarray(sim, float), np.asarray(obs, float)
    mu_s, mu_o = sim.mean(), obs.mean()
    sd_s, sd_o = sim.std(ddof=0), obs.std(ddof=0)
    if mu_o == 0.0 or sd_o == 0.0:
        return float('nan'), float('nan'), float('nan'), float('nan')
    r = float(np.corrcoef(sim, obs)[0, 1])
    alpha = sd_s / sd_o          # variability ratio
    beta = mu_s / mu_o           # bias ratio
    kge_val = 1.0 - np.sqrt((r - 1.0) ** 2 + (alpha - 1.0) ** 2 + (beta - 1.0) ** 2)
    return kge_val, r, alpha, beta


def baseline_metrics(solver_path):
    """Compute NSE and KGE for each baseline solver vs the reference DAE."""
    out = []
    with h5py.File(solver_path, 'r') as h:
        ref_Q = np.array(h['comparison/ref_Q']).squeeze()
        ref_t = np.array(h['comparison/ref_t']).squeeze()
        ref_label = read_str(h, h['comparison/ref_label'])
        n_solv = int(np.array(h['comparison/n_solvers']).squeeze())
        Q_refs = h['comparison/solvers/Q_outlet']
        t_refs = h['comparison/solvers/t']
        lab_refs = h['comparison/solvers/label']
        nse_refs = h['comparison/solvers/nse']
        for i in range(n_solv):
            Q_b = np.array(h[Q_refs[i, 0]]).squeeze()
            t_b = np.array(h[t_refs[i, 0]]).squeeze()
            lab = read_str(h, h[lab_refs[i, 0]])
            nse_stored = float(np.array(h[nse_refs[i, 0]]).squeeze())
            # Match the stored NSE convention used in the manuscript:
            # obs = DAE reference, sim = baseline solver, evaluated on the
            # baseline time grid (verified to reproduce nse_stored to 4 d.p.).
            dae_on_b = np.interp(t_b, ref_t, ref_Q)
            nse_v = nse(Q_b, dae_on_b)
            kge_v, r, alpha, beta = kge(Q_b, dae_on_b)
            out.append(dict(label=lab, nse_stored=nse_stored,
                            nse=nse_v,
                            kge=kge_v, r=r, alpha=alpha, beta=beta))
    return ref_label, out


# --------------------------------------------------------------------------
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
CASES = [
    ('V-Tilted',     'vtilted',      'vtilted'),
    ('Walnut Gulch', 'walnut_gulch', 'walnut_gulch'),
]


def _case_paths(case_dir, file_prefix):
    base = os.path.join(ROOT, 'results', case_dir, 'AWR2026')
    return dict(
        mc=os.path.join(base, f'{file_prefix}_mc_validation_ca_routing_ALL.mat'),
        solv=os.path.join(base, f'{file_prefix}_solver_comparison.mat'),
        out=os.path.join(base, f'{file_prefix}_R1_response_postprocess.mat'),
    )


def compute_c2(mc_path):
    t, Q_mean, Q_std, Q_ens = load_mc_outlet(mc_path)
    ipk = int(np.argmax(Q_mean))
    x_peak = Q_ens[:, ipk]
    mu_pk = float(Q_mean[ipk])
    sd_pk = float(Q_std[ipk])
    g1, g2, cov = skewness_kurtosis_coverage(x_peak, mu_pk, sd_pk, k=1.96)
    return dict(
        peak_index=int(ipk),
        peak_time=float(t[ipk]),
        n_MC=int(x_peak.size),
        Q_peak_ensemble=np.asarray(x_peak, float).reshape(-1, 1),
        mu_peak=mu_pk,
        sigma_peak=sd_pk,
        skewness_gamma1=float(g1),
        excess_kurtosis_gamma2=float(g2),
        empirical_coverage=float(cov),
        nominal_coverage=0.95,
        k_sigma=1.96,
        variable='Q_outlet',
        description='C2-R1: shape statistics of MC outlet discharge at peak '
                    'time and empirical coverage of mu +/- 1.96 sigma interval.',
    )


def compute_c3(solv_path):
    ref_label, rows = baseline_metrics(solv_path)
    n = len(rows)
    labels = np.empty((n,), dtype=object)
    nse_v = np.zeros(n); nse_s = np.zeros(n); kge_v = np.zeros(n)
    r_v = np.zeros(n); a_v = np.zeros(n); b_v = np.zeros(n)
    for i, row in enumerate(rows):
        labels[i] = row['label']
        nse_v[i] = row['nse']; nse_s[i] = row['nse_stored']
        kge_v[i] = row['kge']; r_v[i] = row['r']
        a_v[i] = row['alpha']; b_v[i] = row['beta']
    return dict(
        reference_label=ref_label,
        baseline_labels=labels,
        NSE=nse_v.reshape(-1, 1),
        NSE_stored=nse_s.reshape(-1, 1),
        KGE=kge_v.reshape(-1, 1),
        r=r_v.reshape(-1, 1),
        alpha=a_v.reshape(-1, 1),
        beta=b_v.reshape(-1, 1),
        convention='sim = baseline solver on its native t grid;  '
                   'obs = DAE reference linearly interpolated to that grid.',
        description='C3-R1: NSE and KGE of baseline solvers against DAE '
                    'reference, plus KGE decomposition (r, alpha, beta).',
    )


def save_case(name, paths):
    c2 = compute_c2(paths['mc'])
    c3 = compute_c3(paths['solv'])
    mdict = {
        'c2_r1': c2,
        'c3_r1': c3,
        'meta': dict(
            response='AWR2026 R1 -- C2-R1 and C3-R1',
            case=name,
            sources=dict(
                mc_file=os.path.basename(paths['mc']),
                solver_file=os.path.basename(paths['solv']),
            ),
        ),
    }
    savemat(paths['out'], mdict, do_compression=True, oned_as='column')
    print(f'  wrote {paths["out"]}')


def main():
    print('=' * 72)
    print(' AWR2026 R1 -- C2-R1 (MC shape stats and coverage at outlet peak)')
    print('=' * 72)
    for name, case_dir, file_prefix in CASES:
        paths = _case_paths(case_dir, file_prefix)
        c2 = compute_c2(paths['mc'])
        print(f'\n[{name}]')
        print(f'  peak time index = {c2["peak_index"]}  '
              f'({c2["peak_time"]:.1f} s);  n_MC = {c2["n_MC"]}')
        print(f'  MC mean / std at peak: '
              f'{c2["mu_peak"]:.4g} / {c2["sigma_peak"]:.4g}')
        print(f'  skewness   gamma_1       = {c2["skewness_gamma1"]:+.3f}')
        print(f'  excess kurtosis gamma_2  = {c2["excess_kurtosis_gamma2"]:+.3f}')
        print(f'  empirical coverage of mu +/- 1.96 sigma (nominal 95%) = '
              f'{100.0 * c2["empirical_coverage"]:.1f}%')

    print('\n' + '=' * 72)
    print(' AWR2026 R1 -- C3-R1 (NSE and KGE of DAE vs baseline solvers)')
    print('=' * 72)
    for name, case_dir, file_prefix in CASES:
        paths = _case_paths(case_dir, file_prefix)
        _, rows = baseline_metrics(paths['solv'])
        print(f'\n[{name}]')
        for r in rows:
            print(f'  baseline = {r["label"]:>32s}'
                  f' | NSE = {r["nse"]:+.4f}  (stored {r["nse_stored"]:+.4f})'
                  f' | KGE = {r["kge"]:+.4f}'
                  f' | r = {r["r"]:+.4f},'
                  f' alpha = {r["alpha"]:+.4f},'
                  f' beta = {r["beta"]:+.4f}')

    print('\n' + '=' * 72)
    print(' Saving per-watershed .mat files')
    print('=' * 72)
    for name, case_dir, file_prefix in CASES:
        save_case(name, _case_paths(case_dir, file_prefix))


if __name__ == '__main__':
    main()
