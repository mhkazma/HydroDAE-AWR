# HydroDAE-AWR

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24545573.svg)](https://doi.org/10.5281/zenodo.24545573)

**Companion code and data for:**

> **Exploring Uncertainty Propagation in Coupled Hydrologic and Hydrodynamic Systems via Distribution-Agnostic State Space Analysis**
>
> M. H. Kazma and A. F. Taha
>
> *Advances in Water Resources*, 215, 105389 (2026). Submitted 5 March 2026; accepted 15 June 2026; published 23 June 2026. https://doi.org/10.1016/j.advwatres.2026.105389 (arXiv: https://arxiv.org/abs/2603.05740)

---

## Overview

Flood forecasts depend on uncertain rainfall, soil properties, and initial conditions. The paper writes two-dimensional overland flow routing and Green-Ampt infiltration as a coupled differential algebraic equation (DAE) state space model, in which the water depth and cumulative infiltration of each cell are differential states and the discharges between cells are algebraic variables given by Manning's equation. On this model, it propagates the uncertainty of rainfall and soil parameters and estimates the watershed states from partial measurements. The approach is distribution-agnostic, since it needs only the covariances of the uncertainties and not their distributions, and it gives confidence intervals for water depths and discharges at measured and unmeasured locations. The code reproduces the experiments of the paper on the synthetic V-tilted catchment and on the Walnut Gulch Experimental Watershed in Arizona, including the comparison with Monte Carlo ensembles of cellular automata and local inertial solvers.

This repository contains the MATLAB simulation code and pre-computed results for reproducing the numerical experiments in the paper.

Two test cases are included:

| Case | Description | Grid | Resolution |
|------|-------------|------|------------|
| **V-Tilted** | Synthetic V-shaped catchment | 50 × 81 | 20 m |
| **Walnut Gulch** | USDA-ARS experimental watershed, Arizona | 242 × 521 | 50 m |

## Repository Structure

```
HydroDAE-AWR/
├── run_AWR2026_paper.m              # Main driver script
├── params_vtilted_AWR2026.m         # V-Tilted watershed configuration
├── params_walnut_gulch_AWR2026.m    # Walnut Gulch watershed configuration
│
├── src/
│   ├── core/                        # DAE solvers
│   │   ├── run_manning_ca_bdf.m           # BDF (primary solver)
│   │   ├── run_manning_ca_bdf_newton.m    # BDF 
│   │   ├── run_manning_ca_ode15s.m        # Semi-explicit DAE via ode15s
│   │   └── run_semi_implicit_simulation.m # Semi-implicit local inertial
│   │
│   ├── domain/                      # Domain setup
│   │   ├── setup_vtilted_domain.m         # Generate synthetic V-tilted domain
│   │   ├── setup_walnut_gulch_domain.m    # Load Walnut Gulch from processed DEM
│   │   ├── fill_dem_pits.m                # DEM pit-filling (Planchon-Darboux)
│   │   └── add_infiltration_params.m      # Append Green-Ampt parameters
│   │
│   ├── estimation_uncertainty/      # PSE with measurements
│   │   ├── run_pse_with_measurements.m    # Cell-by-cell PSE (primary)
│   │   ├── run_pse_global_sparse.m        # Global sparse PSE (full Jacobian)
│   │   ├── compute_kxx_measured.m         # State covariance via augmented system
│   │   ├── configure_sensors.m            # Sensor placement strategies
│   │   ├── build_measurement_schedule.m   # Temporal measurement scheduling
│   │   └── build_measurement_matrix.m     # Local C and K_vv construction
│   │
│   ├── uncertainty/                 # Open-loop PSE & Monte Carlo
│   │   ├── run_pse_dae.m                  # PSE without measurements (baseline)
│   │   └── run_monte_carlo_validation.m   # MC validation with LHS sampling
│   │
│   └── validation/                  # Comparison solvers
│       └── run_ca_routing_comparison.m    # Original CA Routing (Gomes Jr. 2023)
│
├── data/
│   └── walnut_gulch/processed/      # Pre-processed Walnut Gulch data
│       ├── walnut_gulch_grid_242x521_burned.mat  # DEM + Manning's n + soils
│       ├── walnut_gulch_spatial_data.mat          # Watershed bounds
│       ├── vegetation_parameters.csv              # 24 vegetation → Manning's n
│       ├── ssurgo_parameters.csv                  # 26 soils → Ks, ψ, θ_s
│       └── instrument_grid_mapping.csv            # Gage → grid cell mapping
│
└── results/                         # Pre-computed simulation outputs
    ├── vtilted/AWR2026/
    │   ├── vtilted_dae_simulation.mat
    │   ├── vtilted_solver_comparison.mat
    │   ├── vtilted_pse_measured_ALL.mat
    │   ├── vtilted_pse_openloop_ALL.mat
    │   ├── vtilted_mc_validation_ca_routing_ALL.mat
    │   └── vtilted_mc_validation_semi_implicit_ALL.mat
    └── walnut_gulch/AWR2026/
        ├── walnut_gulch_dae_simulation.mat
        ├── walnut_gulch_solver_comparison.mat
        ├── walnut_gulch_pse_measured_ALL.mat
        ├── walnut_gulch_pse_openloop_ALL.mat
        └── walnut_gulch_mc_validation_ca_routing_ALL.mat
```

## Requirements

- **MATLAB R2020b** or later
- No additional toolboxes are required for the core solvers
- Parallel Computing Toolbox (optional, for Monte Carlo with `MC_USE_PARALLEL = true`)

## Quick Start

### Run with pre-computed results

The `results/` directory contains all simulation outputs. You can load any `.mat` file directly in MATLAB without running the pipeline:

```matlab
load('results/walnut_gulch/AWR2026/walnut_gulch_dae_simulation.mat');
% Variables: t_out, x_out, sim_results, sim_params, rainfall, time, ...
```

### Run the full pipeline

```matlab
% From the repository root:
run_AWR2026_paper
```

This runs all phases for Walnut Gulch by default. To configure:

```matlab
% Select watershed(s)
WATERSHEDS = {'vtilted'};                  % or {'walnut_gulch'}, or both
FORCE_RERUN = false;                       % Skip phases with existing outputs

% Toggle phases
RUN_DAE_SIMULATION = true;                 % Phase 2a: Primary DAE solver
RUN_SOLVER_COMPARISON = true;              % Phase 2b: CA Routing + semi-implicit
RUN_PSE_WITH_MEASUREMENTS = true;          % Phase 3:  PSE uncertainty propagation
RUN_PSE_OPEN_LOOP = true;                  % Phase 3:  PSE without measurements
RUN_MONTE_CARLO = false;                   % Phase 3b: MC validation (slow)

run_AWR2026_paper
```

## Configuration Flags

| Flag | Default | Description |
|------|---------|-------------|
| `WATERSHEDS` | `{'walnut_gulch'}` | Watershed(s) to process |
| `FORCE_RERUN` | `true` | Overwrite existing `.mat` outputs |
| `RUN_DAE_SIMULATION` | `true` | Run primary BDF DAE solver |
| `RUN_SOLVER_COMPARISON` | `true` | Run CA Routing + semi-implicit for comparison |
| `RUN_PSE_WITH_MEASUREMENTS` | `true` | PSE with sensor measurements |
| `RUN_PSE_OPEN_LOOP` | `true` | PSE without measurements (baseline) |
| `PSE_METHOD` | `'cell_by_cell'` | PSE method: `'cell_by_cell'` or `'global_sparse'` |
| `RUN_MONTE_CARLO` | `false` | Monte Carlo validation |
| `MC_N_SAMPLES` | `100` | Number of MC realizations |
| `MEASUREMENT_INTERVAL` | `60` | Sensor measurement interval (seconds) |

## DAE Formulation

Each grid cell carries 6 state variables: **[h, Q_l, Q_r, Q_u, Q_d, F]**

| State | Type | Description |
|-------|------|-------------|
| h | Differential | Water depth (m) |
| Q_l, Q_r | Algebraic | Left/right discharge (m³/s) via Manning's equation |
| Q_u, Q_d | Algebraic | Up/down discharge (m³/s) via Manning's equation |
| F | Differential | Cumulative infiltration (m), Green-Ampt model |

The system is a semi-explicit Index-1 DAE: differential equations govern mass balance (h) and infiltration (F), while algebraic constraints enforce Manning's discharge relations.

## Output File Contents

Each `*_dae_simulation.mat` contains:

| Variable | Description |
|----------|-------------|
| `t_out` | Time vector (s) |
| `x_out` | Full state trajectory [n_times × 6·n_cells] |
| `sim_results` | Struct: `h_final`, `F_final`, `Q_outlet_ts`, `water_balance`, `computation_time` |
| `sim_params` | Domain parameters, grid topology, soil/hydraulic properties |

Each `*_pse_measured_ALL.mat` contains:

| Variable | Description |
|----------|-------------|
| `pse_meas` | Struct: `h_std`, `h_lower`, `h_upper`, `F_std`, `Q_*_std` [ny × nx × nt] |
| `pse_case` | Uncertainty configuration (CVs for rainfall, Manning's n, Ks, ψ, η) |
| `sensor_config` | Sensor placement and noise parameters |

## Citation

```bibtex
@article{kazma2026uncertainty,
  title   = {Exploring uncertainty propagation in coupled hydrologic and hydrodynamic systems
             via distribution-agnostic state space analysis},
  author  = {Kazma, Mohamad H. and Taha, Ahmad F.},
  journal = {Advances in Water Resources},
  volume  = {215},
  pages   = {105389},
  year    = {2026},
  doi     = {10.1016/j.advwatres.2026.105389},
  note    = {Submitted 5 March 2026; accepted 15 June 2026; published 23 June 2026. arXiv:2603.05740}
}
```

To cite the code itself (version 1.0, archived on Zenodo):

```bibtex
@misc{kazma2026hydrodaeawr,
  author       = {Kazma, Mohamad H. and Taha, Ahmad F.},
  title        = {{HydroDAE-AWR}: Uncertainty propagation in coupled hydrologic and hydrodynamic
                  systems via distribution-agnostic state space analysis (Version 1.0)},
  howpublished = {Zenodo},
  year         = {2026},
  doi          = {10.5281/zenodo.24545573}
}
```

## Data Acknowledgment

The Walnut Gulch Experimental Watershed data used in this study were obtained from the USDA-ARS Southwest Watershed Research Center Data Access Project:

> <https://www.tucson.ars.ag.gov/dap/>

All data available through the SWRC data access website are in the public domain and are not restricted by copyright.

> Datasets were provided by the USDA-ARS Southwest Watershed Research Center. Funding for these datasets was provided by the United States Department of Agriculture, Agricultural Research Service.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
