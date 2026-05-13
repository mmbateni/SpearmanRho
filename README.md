# Trend Analysis & Lag Matrix Tools (R)

## Overview
Collection of R functions for autoregressive lag matrix construction, rank-based autocorrelation estimation, and Spearman's rho trend testing with serial correlation correction.

## Scripts
| File | Function | Purpose |
|------|----------|---------|
| `newlagmatrix.R` | `newlagmatrix(x, nlags, c)` | Constructs lag matrix and response vector for autoregressive modeling. |
| `ac1ranks.R` | `ac1ranks(x)` | Computes lag-1 autocorrelation of data ranks. |
| `spearman_rho.R` | `spearman_rho(v, alpha, ac_correction)` | Standard Spearman's rho trend test. |
| `spearman_rho_modified.R` | `spearman_rho_modified(v, alpha, ac_correction)` | Modified Spearman's rho with DRS variance correction for serially correlated data. |
| `calvar_DRS.R` | `calvar_DRS(n, r)` | Computes DRS variance correction factor. |

## Dependencies
Base R only. Relies on `stats::acf`, `stats::rank`, `stats::pnorm`, `stats::qnorm`.

## Usage
```R
source("spearman_rho_modified.R")
result <- spearman_rho_modified(my_timeseries, alpha = 0.05, ac_correction = TRUE)
print(result$trend)
