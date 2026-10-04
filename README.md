# Optimal Currency Bands

Replication code and thesis text for **"Optimal Currency Bands"**.

**Author:** Iana Sekerina

**Supervisor:** Prof. Francesco Lippi

**Department:** Economics and Finance, LUISS Guido Carli, Rome, 2026

---

## Abstract

This thesis provides an analytical solution to the steady-state optimal currency-band problem, which prior target-zone literature addressed only numerically. Building on Krugman (1991), a monetary authority chooses a symmetric exchange-rate band to minimize steady-state exchange-rate variance plus proportional FX-intervention costs. The thesis derives an exact closed-form expression for the steady-state variance, solves the authority's problem analytically to obtain closed-form optimal bands in narrow- and wide-band limits, and shows that the optimal band depends only on the ratio of intervention cost to fundamental volatility.

The model is calibrated to daily Italian lira data from the ERM I period, spanning the lira's 1990 transition from a ±6% to a ±2.25% band. Structural estimates of fundamental volatility and intervention cost are recovered for each regime from the model's exact equilibrium conditions. A cross-regime counterfactual decomposes the observed narrowing of the band into a volatility channel and a cost channel, with volatility accounting for the larger share — though neither channel alone rationalizes the observed ±2.25% width.

## Repository Contents

| File | Description |
|---|---|
| `OptimalBands_Code.m` | MATLAB replication script. Reproduces all figures and calibration results (Tables 3–4, Figures 1–14). |
| `OptimalCurrencyBands_Oct26.pdf` | Full thesis text (theory, proofs, empirical application, appendix). |
| `ert_h_eur_d__custom_19839545_linear.csv` | ECB daily ITL/ECU exchange rate data (main file, not included — see Data below). |
| `ert_h_eur_d__custom_19839614_linear.csv` | ECB daily ITL/ECU exchange rate data (supplementary tail, not included — see Data below). |

## Requirements

- MATLAB R2019b or later (uses `tiledlayout`, `xline`, `yline`)
- No additional toolboxes required — all root-finding uses built-in `fzero`

## Data

The two required ECB Statistical Data Warehouse CSV exports (daily ITL/ECU exchange rate, columns G:H) are not redistributed in this repository. To reproduce the results:

1. Download the daily ITL/ECU series from the [ECB Statistical Data Warehouse](https://sdw.ecb.europa.eu/).
2. Place both CSV files in a single folder.
3. Update `DATA_DIR` near the top of `OptimalBands.m` to point to that folder.

## How to Run

```matlab
% 1. Set DATA_DIR in OptimalBands.m to the folder containing the two ECB CSV files
% 2. Run the script
OptimalBands
```

All figures are generated automatically; calibration results print to the console and are also saved as structured output variables in the workspace:

| Variable | Contents |
|---|---|
| `res_6pct` | Calibration results for the ±6% ERM band (1987–1990) |
| `res_2pct` | Calibration results for the ±2.25% ERM band (1990–1992) |
| `cf1`, `cf2` | Counterfactual forward-solve results (cross-regime decomposition) |

## Key Results (October 2026)

| | ±6% regime | ±2.25% regime |
|---|---|---|
| Estimation window | 12 Jan 1987 – 4 Jan 1990 | 5 Jan 1990 – 11 Sep 1992 |
| Observed band | 6.00% | 2.25% |
| Estimated 𝜎̂ | 0.0612 | 0.0362 |
| Estimated 𝜅̂ | 0.0289 | 0.0171 |
| 𝜃 = 𝜅̂/𝜎̂ | 0.472 | 0.472 |

## Acknowledgments

AI tools were used solely for code debugging and editorial/LaTeX assistance. All research questions, theoretical derivations, empirical analyses, and conclusions are the author's own.

