# Code Review — `asym_stdev_v4_note.pdf` Pipeline

**Date:** 2026-05-07  
**Scope:** All scripts that generate outputs included in `asym_stdev_v4_note.tex` / `.pdf`

---

## Scripts Reviewed

| Script | Role |
|---|---|
| `2a_regression_fm_asym_v4.R` | Static FM regression; produces `fm_stdev_posneg_rebuilt_v4.RDS` |
| `2a_regression_fm_asym_v4_dynamic.R` | Dynamic FM regression; produces `fm_dynamic_stdev_posneg_rebuilt_v4_unbalanced.RDS` |
| `test_v4_delta_posneg.R` | Positive-vs-negative delta tests; produces `v4_delta_posneg_tests_static/dynamic.{RDS,dta,csv}` |
| `make_asym_stdev_v4_note_table.R` | Formats static regression table → `reg_asym_stdev_posneg_rebuilt_v4.tex` |
| `make_asym_stdev_v4_dynamic_note_table.R` | Formats dynamic regression table → `reg_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex` |
| `make_asym_stdev_v4_dynamic_summary_table.R` | Formats dynamic bin summary → `summary_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex` |
| `make_v4_delta_posneg_test_tables.R` | Formats delta test tables → `test_v4_delta_posneg_{static,dynamic}_spec3.tex` |
| `pres_plots_static_asym_pooled.py` | Static presentation plots → `plots/v4/asym_v4_{BMI,FIT,OFI}.png` |
| `pres_plots_dynamic_asym_pooled.py` | Dynamic presentation plots → `plots_dynamic/v4_unbalanced/dynamic_asym_v4_{FIT,OFI}_4lag.png` |

---

## Findings

### Bug — affects generated output

#### 1. Mislabeled column header in both delta test tables

**File:** `make_v4_delta_posneg_test_tables.R`, lines 58 and 74  
**Outputs affected:** `test_v4_delta_posneg_static_spec3.tex`, `test_v4_delta_posneg_dynamic_spec3.tex`

The first column of both tables is labeled `"Shock"` in the header, but the rows contain demand types (`BMI`, `FIT`, `OFI`), not shock signs.

```r
# Line 58 — static table header (wrong)
"Shock & Contrast & $\Delta(+)$ & $\Delta(-)$ & $\Delta(+) - \Delta(-)$ & SE & $t$ & Months \\"

# Line 74 — dynamic table header (wrong)
"Shock & Lag & Contrast & $\Delta(+)$ & $\Delta(-)$ & $\Delta(+) - \Delta(-)$ & SE & $t$ & Months \\"
```

The row functions emit `r$type` / `r$demand_type` as the first cell, which take values `"BMI"`, `"FIT"`, `"OFI"`. The correct header for both tables is `"Type"`.

**Fix:** Replace `"Shock"` with `"Type"` on lines 58 and 74, then re-run `make_v4_delta_posneg_test_tables.R` and recompile the PDF.

---

### Flag — verify against main paper convention

#### 2. Standard FM SE used for stars, not Newey-West SE

**Files:** `make_asym_stdev_v4_note_table.R`, `make_asym_stdev_v4_dynamic_note_table.R`

Both table generators compute significance stars using `r$se` (standard FM standard error):

```r
fmt_coef <- function(coef, se) {
  n_stars <- sum(abs(coef / se) > cuts)
  ...
}
# called as: fmt_coef(r$coef, r$se)
```

The FM helper function (`p.fama_macbeth_signbins`) computes both `se` (standard: `sd/sqrt(T)`) and `se_nw` (Newey-West corrected), but only `se` is passed to the table formatter. The PLAN says star thresholds should match `reg_static_fm_stdev.Rmd` — **check whether that script uses `se` or `se_nw` for stars**. If the main paper tables use NW-corrected SE, the v4 tables are inconsistent.

---

### Minor fragility — works correctly but worth knowing

#### 3. Positional column assignment in dynamic summary table

**File:** `make_asym_stdev_v4_dynamic_summary_table.R`, lines 46–50

```r
summary_data <- data_reg[, .(
  shock = fcase(cumofi_lag > 0, "Positive", cumofi_lag < 0, "Negative", default = NA_character_),
  bin   = fcase(...)
)]
summary_data[, `:=`(
  type       = data_reg$type,
  cumofi_lag = data_reg$cumofi_lag,
  ofi        = data_reg$ofi
)]
```

`summary_data` is created with only `shock` and `bin`, then `type`, `cumofi_lag`, `ofi` are attached by positional alignment from `data_reg`. This is correct because `data_reg[, .(col1, col2)]` preserves row order in data.table. But it would silently break if rows were reordered between the two statements. The safer pattern is to compute all columns in one call:

```r
summary_data <- data_reg[, .(
  shock      = fcase(...),
  bin        = fcase(...),
  type       = type,
  cumofi_lag = cumofi_lag,
  ofi        = ofi
)]
```

---

## Verification Status

| Check | Method | Status |
|---|---|---|
| Bin boundary logic matches LaTeX spec | Static code read | ✅ Verified |
| Bins mutually exclusive, cover all non-zero obs | Static code read | ✅ Verified |
| Difference rows use common-month inner join | Static code read | ✅ Verified |
| FM SE formula consistent across scripts | Static code read | ✅ Verified |
| Star thresholds match PLAN conventions | Static code read | ✅ Verified |
| `var_added` filter selects spec_idx=3 for plots/tables | Static code read | ✅ Verified |
| Dynamic bins use `cumofi_lag` for membership, `ofi` as regressor | Static code read | ✅ Verified |
| SD grouping equivalent between regression and summary scripts | Static code read | ✅ Verified |
| v4 dynamic script mirrors reference dynamic script architecture | Read reference + v4 side-by-side | ✅ Verified |
| Note text highlights match actual RDS outputs | Loaded RDS, checked all cited values | ✅ Verified |
| Table 1 bin summary statistics match actual data | Recomputed from raw data | ✅ Verified |

---

## Verified Correct

The following were checked thoroughly and are implemented correctly.

### Bin construction

**Static (`2a_regression_fm_asym_v4.R` lines 201–207):**

```r
ofi_bin1_pos := ofi * (ofi > 0 & ofi < ofi_sd_pos)       # 0 < d < σ+
ofi_bin2_pos := ofi * (ofi >= ofi_sd_pos & ofi < 2*ofi_sd_pos)  # σ+ ≤ d < 2σ+
ofi_bin3_pos := ofi * (ofi >= 2 * ofi_sd_pos)             # d ≥ 2σ+

ofi_bin1_neg := ofi * (ofi < 0 & ofi > -ofi_sd_neg)       # -σ- < d < 0
ofi_bin2_neg := ofi * (ofi <= -ofi_sd_neg & ofi > -2*ofi_sd_neg) # -2σ- < d ≤ -σ-
ofi_bin3_neg := ofi * (ofi <= -2 * ofi_sd_neg)            # d ≤ -2σ-
```

- Boundary at ±σ belongs to bin 2, not bin 1. Matches the LaTeX definition exactly.
- Bins are mutually exclusive and cover all non-zero observations.
- Zero-demand rows (`ofi == 0`) contribute zero to all six regressors (they land in no bin and affect only the intercept). Appropriate.
- Sign-specific SDs are computed within `(type, yyyymm)`. Validated: zero missing SD months for all three demand types.

**Dynamic (`2a_regression_fm_asym_v4_dynamic.R` lines 142–148):**

- Bin membership correctly based on `cumofi_lag` (lagged cumulative demand), not contemporaneous `ofi`.
- Regressor value inside each bin is contemporaneous `ofi`. Matches the note's dynamic specification.
- SD grouping uses `by = .(yyyymm, type)` where `type` encodes both demand type and horizon (`"FIT_4lag"` etc.) — equivalent to the summary script's `by = .(yyyymm, type, hor)`. ✓

### Difference rows use the common-month sample

In `p.fama_macbeth_signbins`, pairwise differences are computed by merging monthly bin coefficients with `all = FALSE` — months where either coefficient is missing are excluded from that difference's FM average. This is the correct approach and is documented in the note. The same logic holds in `compute_delta_tests` (`test_v4_delta_posneg.R`): `complete.cases` on the four needed bin columns. ✓

### FM standard error formula is consistent across scripts

`p.fama_macbeth_signbins` uses `sqrt(vcov(lm(coef~1))[1,1])` = `sd/sqrt(T)`. The `fm_mean` helper in `test_v4_delta_posneg.R` uses `sd(x)/sqrt(length(x))`. Identical formula. ✓

### Significance star thresholds

`cuts <- abs(qnorm(c(.01, .05, .1) / 2))` → `[2.576, 1.960, 1.645]`. `sum(|t| > cuts)` gives 0/1/2/3 stars at 10%/5%/1% two-sided significance. Matches PLAN conventions. ✓

### `var_added` filter for most-controlled spec

Both plot scripts filter `df["var_added"] == "controls_char+controls_liq"` to select the most-controlled specification. The `spec_labels` table in both regression scripts maps `spec_idx == 3` to exactly this label. ✓

### Progressive interaction specs carry the correct base formula

The second loop in `p.process_one_type` appends interaction terms to `ff` from `spec_idx = 3` (the full char+liq controls spec, plus BMI-specific controls for BMI). The R for-loop leaves `spec_idx = 3` after the first loop; the second loop increments from 4 onward. This matches the v2 architecture. ✓

### Table column-count declarations match output

| Table | Declaration | Cols |
|---|---|---|
| Static regression | `{lccccccccc}` | 10 (label + 9 data) |
| Dynamic regression | `{lcccccc}` | 7 (label + 6 data) |
| Delta test static | `{llrrrrrr}` | 8 |
| Delta test dynamic | `{llrrrrrrr}` | 9 |

All match their row-building functions. ✓

### `\vspace` in generated tables

`\vspace{5pt}` lines appear as standalone strings inside the tabular body. The note wraps these inputs in `\begingroup \renewcommand{\vspace}[1]{} ... \endgroup`, which suppresses them harmlessly. No LaTeX errors. ✓

### Dynamic regression test script re-runs bins consistently

`test_v4_delta_posneg.R` rebuilds bins from scratch for both static and dynamic tests, using the same RDS inputs and the same SD grouping as the regression scripts. The monthly coefficients from the test script match those from the regression script for `spec_idx = 1, 2, 3`. ✓

---

## Action Items

| Priority | Action | File |
|---|---|---|
| **Fix** | Replace `"Shock"` → `"Type"` on lines 58 and 74; re-run and recompile | `make_v4_delta_posneg_test_tables.R` |
| **Verify** | Check whether `reg_static_fm_stdev.Rmd` uses `se` or `se_nw` for stars; align v4 tables if needed | `make_asym_stdev_v4_note_table.R`, `make_asym_stdev_v4_dynamic_note_table.R` |
| *Optional* | Consolidate positional assignment into single `[, .()]` call | `make_asym_stdev_v4_dynamic_summary_table.R` L46–50 |
