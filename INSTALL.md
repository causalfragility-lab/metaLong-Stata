# metaLong for Stata 14.1

**Longitudinal Meta-Analysis with Robust Variance Estimation and Sensitivity Analysis**

Stata translation of the R package `metaLong` v0.1.0  
Author: Subir Hait · Michigan State University · haitsubi@msu.edu · ORCID 0009-0004-9871-9677

---

## Installation

### Option A — SSC (recommended)

Once the package is listed on SSC, install with a single command:

```stata
. ssc install metalong, replace
```

This installs all ado-files and help files automatically into your **PLUS**
ado-path, which is the correct location for user-written packages.  See
`help sysdir` for the location of PLUS on your system.

### Option B — net install (from a web server or GitHub)

```stata
. net install metalong, from("https://raw.githubusercontent.com/causalfragility-lab/metalong/main") replace
```

### Option C — Manual install from this ZIP

1. Unzip the archive. You will find ado-files (`.ado`) and help files (`.sthlp`).

2. Find your **PLUS** directory:

   ```stata
   . sysdir
   ```

   Typical PLUS paths:
   - **Windows**: `C:\Users\<you>\ado\plus\`
   - **Mac/Linux**: `~/ado/plus/`

   > **Note:** Install into **PLUS**, not PERSONAL. PERSONAL is intended for
   > your own private programs; PLUS is the correct location for installed
   > packages. See the Statalist discussion by Christopher Baum at
   > https://www.statalist.org/forums/forum/general-stata-discussion/general/1785603

3. Copy all `.ado` files and `.sthlp` files into the appropriate
   subdirectory of PLUS (e.g., `plus/m/` for commands starting with "m").

4. Verify in Stata:

   ```stata
   . which ml_meta
   . help metalong
   ```

---

## Quick-start workflow

```stata
* Set your working directory first
cd "path/to/your/project"

* Step 1: Simulate longitudinal meta-analytic data
sim_longmeta, k(20) times(0 6 12 24) seed(42) clear

* Step 2: Pool effects at each time point
ml_meta yi vi, study(study) time(time) saving(meta_res) replace

* Step 3: ITCV sensitivity analysis
ml_sens yi vi, study(study) time(time) ///
    metafile(meta_res) delta(0.15) saving(sens_res) replace

* Step 4: Benchmark ITCV against observed covariates
ml_benchmark yi vi, study(study) time(time) ///
    metafile(meta_res) sensfile(sens_res) ///
    covariates(pub_year quality n) saving(bench_res) replace

* Step 5: Leave-k-out fragility analysis
ml_fragility yi vi, study(study) time(time) ///
    metafile(meta_res) maxk(3) saving(frag_res) replace

* Step 6: Spline time trend
ml_spline, metafile(meta_res) df(3) saving(spline_res) replace

* Step 7: Combined figure
metalong_plot, metafile(meta_res) sensfile(sens_res)    ///
    splinefile(spline_res) fragfile(frag_res)            ///
    saving(figure.gph) replace

graph export figure.pdf, replace
```

A complete self-contained example is provided in `run_metalong.do`. Set your
working directory before running it — all output files are saved there.

---

## Command reference

| Command          | Purpose                                               |
|------------------|-------------------------------------------------------|
| `sim_longmeta`   | Simulate a longitudinal meta-analytic dataset         |
| `ml_meta`        | Pool effects at each time point (DL + RVE)           |
| `ml_sens`        | Time-varying ITCV sensitivity analysis                |
| `ml_benchmark`   | Benchmark ITCV against observed covariates            |
| `ml_spline`      | Restricted cubic spline time trend                    |
| `ml_fragility`   | Leave-k-out fragility index                           |
| `metalong_plot`  | Combined publication figure                           |

Type `help metalong` in Stata for the full package overview.

---

## Statistical methods summary

### Pooling (`ml_meta`)
Uses the **DerSimonian-Laird** (DL) estimator of between-study variance τ².  
Random-effects weights: *wᵢ* = 1/(vᵢ + τ²).  
Standard errors are **cluster-robust** (Huber-White sandwich, clustered by study),  
analogous to the Hedges-Tipton-Johnson (2010) RVE.  
Small-sample correction: **t(k − 1)** degrees of freedom (CR1 approximation).

> **Note on CR2/Satterthwaite**: Stata 14.1 does not have a built-in CR2
> sandwich estimator with Satterthwaite degrees of freedom (Tipton 2015).
> The t(k − 1) approach is the best available native approximation.

### Sensitivity (`ml_sens`)
- sy² = Σwᵢ(yᵢ − θ̂)² / Σwᵢ
- r = θ̂ / √(θ̂² + sy²);  ITCV = √|r|
- θ* = |θ̂| − t_crit · SE;  ITCV_adj = √|r*| (or 0 if θ* ≤ 0)
- Fragile if ITCV_adj < δ (default δ = 0.15)

### Benchmark (`ml_benchmark`)
For each covariate Z at time t:  
r_partial = t_Z / √(t_Z² + df), compared to ITCV_adj(t).

### Fragility (`ml_fragility`)
LOO + random leave-*k*-out search (up to 500 combinations per *k*).

---

## Requirements

- **Stata 14.1** or later (no additional packages required)
- `mkspline` (built into Stata) used by `ml_spline`
- `graph combine` (built into Stata) used by `metalong_plot`

---

## References

- Frank, K. A. (2000). Impact of a confounding variable on a regression
  coefficient. *Sociological Methods & Research*, 29(2), 147–194.
  doi:10.1177/0049124100029002003

- Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
  estimation in meta-regression with dependent effect size estimates.
  *Research Synthesis Methods*, 1(1), 39–65. doi:10.1002/jrsm.5

- Tipton, E. (2015). Small sample adjustments for robust variance
  estimation with meta-regression. *Psychological Methods*, 20(3),
  375–393. doi:10.1037/met0000011

---

## License

MIT License — see LICENSE file.
