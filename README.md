# FO_LEs_ISR

Official MATLAB code release for the paper:

**Unveiling and resolving algorithmic false chaos in fractional-order systems via infinite state representation**

Zhimeng Dong and Yong Xie, *Nonlinear Dynamics* (2026).

## Scope

This compact release contains the fractional Lorenz example at `alpha = 0.995` used in Fig. 4 of the paper. It computes the augmented Lyapunov spectrum with the infinite state representation (ISR), isolates the physical Lyapunov exponents with the degeneracy-rejection algorithm, and performs two independent checks using two-trajectory tracking and the 0-1 test for chaos.

The supplied precomputed data are the data used for the reported figure. With the threshold `delta = 0.1`, the expected results are approximately:

- Physical Lyapunov exponents: `(0.826, -0.017, -14.834)`. The near-zero exponent is reported as `0` in the figure caption.
- Two-trajectory LLE estimate: `0.895`.
- Median 0-1 test statistic: `0.997`.

## Numerical Settings

| Quantity | Value |
|---|---:|
| Fractional order `alpha` | `0.995` |
| Lorenz parameters `(sigma, rho, beta)` | `(10, 28, 8/3)` |
| Frequency nodes `N` | `10` |
| Approximation band | `[1e-4, 1e4]` |
| Transient duration | `50` |
| LE computation duration | `200` |
| QR interval | `0.5` |
| ODE15s relative/absolute tolerances | `1e-10 / 1e-10` |
| Degeneracy threshold `delta` | `0.1` |

## Requirements

- MATLAB R2024a or later is recommended.
- The default Oustaloup calculation uses standard MATLAB functions.
- `generate_optimization_params.m` is included for the alternative kernel approximation described in the paper and requires Optimization Toolbox when used.

## Quick Start

```bash
git clone https://github.com/ZhimengDong/FO_LEs_ISR.git
```

Open MATLAB, change to the repository root, and run:

```matlab
main
```

By default, `main.m` loads the paper data from `Data/Lorenz_ISR_Results_alpha_0.995.mat`, prints the extracted values, and recreates `Lorenz_ISR_Results.png`.

To perform the full computation, run:

```matlab
main(true)
```

The recomputed result is saved under a separate filename, so the supplied paper data are not overwritten. Fractional finite-time Lyapunov estimates may show small numerical differences when recomputed, while the expected chaotic signature remains `(+, 0, -)`.

## Repository Layout

```text
.
|-- main.m
|-- Lorenz_ISR_Results.png
|-- LICENSE
|-- Data/
|   `-- Lorenz_ISR_Results_alpha_0.995.mat
|-- Src/
|   |-- generate_ISR_timeseries.m
|   |-- generate_optimization_params.m
|   |-- generate_oustaloup_params.m
|   `-- run_ISR_LEs.m
`-- Utils/
    |-- compute_01_test.m
    |-- compute_LLE_experimental.m
    |-- extract_physical_LEs.m
    `-- plot_Lorenz_ISR_results.m
```

## Citation

Please cite the published paper when using this code or the supplied data.
