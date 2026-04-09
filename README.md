# metaLong for Stata 14.1

**Longitudinal Meta-Analysis with Robust Variance Estimation and Sensitivity Analysis**

[![Stata](https://img.shields.io/badge/Stata-14.1%2B-blue)](https://www.stata.com)
[![License: MIT](https://img.shields.io/badge/License-MIT-green)](LICENSE)
[![SSC](https://img.shields.io/badge/SSC-metalong-orange)](https://ideas.repec.org/s/boc/bocode.html)

Stata translation of the R package [`metaLong`](https://github.com/causalfragility-lab/metaLong) v0.1.0.  
Author: **Subir Hait** · Michigan State University · haitsubi@msu.edu · [ORCID 0009-0004-9871-9677](https://orcid.org/0009-0004-9871-9677)

---

## Overview

`metaLong` provides a complete workflow for synthesising evidence from studies that report outcomes at **multiple follow-up time points** (longitudinal meta-analysis). It is the only Stata package designed specifically for this setting.

Most meta-analysis packages in Stata treat each time point independently or require manual looping. `metaLong` pools effects across time jointly, propagates uncertainty through sensitivity and robustness analyses, and produces publication-ready figures — all within a single coherent pipeline.

---

## Commands

| Command | Purpose |
|---|---|
| `sim_longmeta` | Simulate a longitudinal meta-analytic dataset |
| `ml_meta` | Pool effects at each time point (DerSimonian-Laird + cluster-robust RVE) |
| `ml_sens` | Time-varying ITCV sensitivity analysis |
| `ml_benchmark` | Benchmark ITCV against observed study-level covariates |
| `ml_spline` | Restricted cubic spline time trend |
| `ml_fragility` | Leave-one-out and leave-k-out fragility index |
| `metalong_plot` | Combined publication-ready figure |

---

## Installation

### From SSC (recommended once listed)

```stata
ssc install metalong, replace
```

### From this repository

```stata
net install metalong, ///
    from("https://raw.githubusercontent.com/causalfragility-lab/metaLong-Stata/main") ///
    replace
```

### Manual install

1. Download or clone this repository
2. Copy all `.ado` and `.sthlp` files to your Stata **PLUS** directory:
   ```stata
   sysdir          /* shows your PLUS path */
   ```
3. Verify:
   ```stata
   which ml_meta
   help metalong
   ```

> **Note:** Install into **PLUS**, not PERSONAL.  
> See [Kit Baum's Statalist post](https://www.statalist.org/forums/forum/general-stata-discussion/general/1785603-use-of-personal-vs-plus) for explanation.

---

## Quick Start

```stata
* Set your working directory
cd "path/to/your/project"

* Step 1: Simulate longitudinal meta-analytic data
sim_longmeta, k(20) times(0 6 12 24) seed(42) clear

* Step 2: Pool effects at each time point
ml_meta yi vi, study(study) time(time) saving(meta_res) replace

* Step 3: ITCV sensitivity analysis
ml_sens yi vi, study(study) time(time) ///
    metafile(meta_res) delta(0.15) saving(sens_res) replace

* Step 4: Benchmark ITCV against observed covariates
ml_benchmark yi vi, study(study) time(time)       ///
    metafile(meta_res) sensfile(sens_res)          ///
    covariates(pub_year quality n) saving(bench_res) replace

* Step 5: Leave-k-out fragility analysis
ml_fragility yi vi, study(study) time(time) ///
    metafile(meta_res) maxk(3) saving(frag_res) replace

* Step 6: Spline time trend
ml_spline, metafile(meta_res) df(3) saving(spline_res) replace

* Step 7: Combined figure
metalong_plot, metafile(meta_res) sensfile(sens_res)  ///
    splinefile(spline_res) fragfile(frag_res)          ///
    saving(figure.gph) replace

graph export figure.pdf, replace
```

A complete self-contained example is in [`run_metalong.do`](run_metalong.do).

---

## Statistical Methods

### Pooling — `ml_meta`

Uses the **DerSimonian-Laird** (DL) estimator of between-study variance τ².  
Random-effects weights: *wᵢ* = 1/(vᵢ + τ²).  
Standard errors are **cluster-robust** (Huber-White sandwich clustered by study),  
analogous to the Hedges, Tipton & Johnson (2010) RVE framework.  
Small-sample correction: **t(k − 1)** degrees of freedom (CR1 approximation).

> **Note on CR2/Satterthwaite:** Stata 14.1 does not have a native CR2 sandwich estimator with Satterthwaite degrees of freedom (Tipton, 2015). The t(k − 1) approach is the best available native approximation.

### Sensitivity — `ml_sens`

At each time point *t*, computes the **Impact Threshold for a Confounding Variable (ITCV)** (Frank, 2000):

- sy² = Σwᵢ(yᵢ − θ̂)² / Σwᵢ
- r = θ̂ / √(θ̂² + sy²) → ITCV = √|r|
- θ\* = |θ̂| − t_crit · SE → ITCV_adj = √|r\*| (or 0 if θ\* ≤ 0)
- A time point is **fragile** if ITCV_adj < δ (default δ = 0.15)

### Benchmark — `ml_benchmark`

For each observed covariate Z at time t, fits a WLS meta-regression and computes:

r_partial = t_Z / √(t_Z² + df)

Compares r_partial against ITCV_adj(t). If r_partial ≥ ITCV_adj, that covariate alone
would be sufficient to nullify the effect.

### Fragility — `ml_fragility`

Implements leave-one-out (LOO) and random leave-*k*-out search
(up to 500 combinations per *k*). Reports the minimum number of study
removals required to flip the significance conclusion at each time point.

### Spline — `ml_spline`

Fits a **restricted cubic spline** meta-regression of pooled effects over
follow-up time using `mkspline`. Produces smooth predicted trajectories
with pointwise confidence bands.

---

## Requirements

- **Stata 14.1** or later
- No additional user-written packages required
- `mkspline` and `graph combine` are both built into Stata

---

## Related Packages

| Package | Language | Description |
|---|---|---|
| [`metaLong`](https://github.com/causalfragility-lab/metaLong) | R (CRAN) | Original R implementation |
| [`drmeta`](https://github.com/causalfragility-lab/drmeta) | R (CRAN) | Design-Robust Meta-Analysis |

---

## References

Frank, K. A. (2000). Impact of a confounding variable on a regression
coefficient. *Sociological Methods & Research*, 29(2), 147–194.
https://doi.org/10.1177/0049124100029002003

Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
estimation in meta-regression with dependent effect size estimates.
*Research Synthesis Methods*, 1(1), 39–65.
https://doi.org/10.1002/jrsm.5

Tipton, E. (2015). Small sample adjustments for robust variance
estimation with meta-regression. *Psychological Methods*, 20(3), 375–393.
https://doi.org/10.1037/met0000011

---

## Citation

If you use `metaLong` for Stata in your research, please cite:

> Hait, S. (2026). *metaLong for Stata 14.1: Longitudinal meta-analysis
> with robust variance estimation and sensitivity analysis.*
> Statistical Software Components, Boston College Department of Economics.
> https://github.com/causalfragility-lab/metaLong-Stata

---

## License

MIT License — see [LICENSE](LICENSE) file.
