# Asymmetry Table — Results

> Mirrors PLAN.md structure. Updated after each step with key findings.
> New agents: read PLAN.md for what to do, RESULTS.md for what was found.

**Last updated:** 2026-04-28 (v4 static/dynamic PDF note refresh)
**Status:** In Progress

---

## Task 1: Generate asymmetry table

**Status:** IMPLEMENTED (awaiting review)

### Key Findings
- Both RDS files have 3,058 rows each; 4 types present (`BMI`, `FIT`, `OFI`, `OFI_pre_whitened`) — `OFI_pre_whitened` correctly excluded from table.
- After filtering to spec_idx ∈ {1,2,3} and difference vars: 36 rows each (9 extra = OFI_pre_whitened, silently ignored in table construction).
- All 27 required combos (3 types × 3 specs × 3 diffs) present with no NaN coef/se.
- OFI spec3 differences (most controlled, pos shocks): bin2−bin1 = −3.30 (se 0.19, ***), bin3−bin1 = −5.67 (se 0.27, ***), bin3−bin2 = −2.37 (se 0.14, ***). Consistent with asymmetry plots showing large negative differences for OFI.
- All differences are negative (multipliers decline with shock size) for both pos and neg shocks, matching prior plots.

### Row Counts / Sample
- Input pos: 3,058 rows; after filter: 36 rows (27 used)
- Input neg: 3,058 rows; after filter: 36 rows (27 used)

### Notes
- The v1 results are from the original asymmetric spec (`2a_regression_fm_asym_AC.R`), which computes σ on zeroed-out OFI (negatives set to 0 but rows kept). This yields a much smaller σ than the pooled spec, making bin assignments non-comparable and multipliers inflated (~2–3× vs pooled). OFI spec3 bin1 multiplier: v1_pos = 8.25 vs pooled = 3.13. Differences are still directionally correct.

### Output
- `code/g_asym/archive/reg_asym_stdev.tex` — 2,472 chars, two-panel LaTeX tabular (differences only; v1 data)

---

## Task 2: Re-run asymmetric FM with pooled σ and row-dropping (v2)

**Status:** APPROVED in direct-mode self-review

### Key Findings
- Script `code/g_asym/2a_regression_fm_asym_v2.R` run successfully.
- **Filtering verified correct** (see validation section below).
- **Positive panel** shows clean declining multipliers for all three types: OFI spec3 bin1/2/3 = 2.99 / 2.27 / 1.56 (all ***). FIT spec3 = 6.43 / 5.13 / 3.57 (all ***). BMI mostly *.
- **Negative panel — OFI**: clean declining pattern, bin3−bin1 = −1.34 (spec3, ***). Consistent with positive panel.
- **Negative panel — FIT**: multipliers much smaller than positive panel (bin1 ≈ 2.47 vs 6.43). Non-monotone: bin3 slightly larger than bin2 (bin3−bin2 = +0.27 to +0.39, all insignificant). Not a bug — see validation section.
- **Negative panel — BMI**: 4,681 obs only; all coefficients insignificant (t << 1.645). Non-monotone pattern is pure noise at this sample size.
- v2 positive multipliers are in line with pooled spec (OFI spec1 bin1 = 2.44 vs pooled 3.13), as expected from using the same pooled σ.

### Filtering Validation (2026-04-26)
Investigated why negative panel results looked strange. Confirmed:
1. **Bin construction** (`1_construct_regression_table_static.R` lines 41–46): `bin := 1` initialized for all rows before overrides, then `ofi_bin{k} = ifelse(bin==k, ofi, 0)`. For negative-OFI rows, the bin variables carry the signed (negative) OFI value. Correct.
2. **`ofi` is type-specific** (`2_put_together_main_data.R` line 78): long-format stack where `ofi` holds BMI demand for BMI rows, FIT demand for FIT rows, OFI demand for OFI rows. The filter `data_all[ofi < 0]` correctly selects negative demand for each type.
3. **Pooled σ is also type-specific**: `sd(ofi)` is computed `by = .(yyyymm, type)`, so BMI bins use BMI σ, FIT bins use FIT σ, etc.
4. **Conclusion**: filtering and bin assignment are both correct. Strange negative results are genuine, driven by small BMI sample and weaker / asymmetric FIT price impact on the sell side.

### Row Counts / Sample
- Positive panel: BMI 5,019 obs, FIT 307,873 obs, OFI 220,873 obs
- Negative panel: BMI 4,681 obs, FIT 234,705 obs, OFI 306,871 obs

### Output
- `code/g_asym/reg_asym_stdev_v2.tex` — generated; two-panel LaTeX table (Panel A: Positive, Panel B: Negative), each with raw multipliers + control rows + differences block

---

## Task 3: Sign-specific σ bins (v3)

**Status:** IMPLEMENTED and run; outputs validated

### Key Findings
- Script `code/g_asym/archive/2a_regression_fm_asym_v3.R` run successfully before archival.
- v3 uses the same row-dropping sign-specific samples as v2, then recomputes `ofi_bin1/2/3` using sign-specific σ by `yyyymm` and `type`.
- Positive panel remains cleanly declining for OFI and FIT. Spec3 OFI bin1/2/3 = 3.23 / 2.63 / 1.72; differences are bin2−bin1 = −0.60 (t = −5.06), bin3−bin2 = −0.91 (t = −9.17), bin3−bin1 = −1.51 (t = −10.76). FIT spec3 bin1/2/3 = 6.54 / 5.70 / 3.96; bin3−bin1 = −2.59 (t = −2.36). BMI remains noisy.
- Negative panel OFI is strongly declining. Spec3 OFI bin1/2/3 = 1.84 / 1.34 / 0.49; differences are bin2−bin1 = −0.51 (t = −6.76), bin3−bin2 = −0.85 (t = −12.68), bin3−bin1 = −1.35 (t = −15.29).
- Negative panel FIT remains non-monotone under sign-specific σ. Spec3 FIT bin1/2/3 = 0.80 / 2.15 / 2.49; bin3−bin1 = +1.69 but insignificant (t = 1.09). The sign-specific binning does not resolve the weak/asymmetric sell-side FIT pattern identified in v2.
- Negative panel BMI remains noisy and insignificant, consistent with the small sample.

### Filtering / Bin Validation
- Row counts after sign filter match v2: positive panel BMI 5,019, FIT 307,873, OFI 220,873; negative panel BMI 4,681, FIT 234,705, OFI 306,871.
- Reconstructed v2 and v3 bin assignments directly from `reg_table_static.RDS`; in every shock/type cell, mean absolute `ofi` rises monotonically from bin 1 to bin 3. The negative FIT non-monotonic multiplier pattern is therefore not caused by reversed or malformed bins.
- v2 pooled-σ negative FIT bin summary: bin 1 = 176,791 obs, mean `ofi` −0.001351, sd 0.001312; bin 2 = 40,711 obs, mean `ofi` −0.005666, sd 0.002131; bin 3 = 17,203 obs, mean `ofi` −0.011023, sd 0.004662.
- v3 sign-specific-σ negative FIT bin summary: bin 1 = 155,449 obs, mean `ofi` −0.001062, sd 0.001012; bin 2 = 49,096 obs, mean `ofi` −0.004405, sd 0.001770; bin 3 = 30,160 obs, mean `ofi` −0.009209, sd 0.004401.
- Full v3 recomputed bin counts by table type:
  - Positive BMI: 3,765 / 670 / 584 across bins 1/2/3.
  - Positive FIT: 223,728 / 52,235 / 31,910 across bins 1/2/3.
  - Positive OFI: 163,014 / 37,776 / 20,083 across bins 1/2/3.
  - Negative BMI: 3,671 / 477 / 533 across bins 1/2/3.
  - Negative FIT: 155,449 / 49,096 / 30,160 across bins 1/2/3.
  - Negative OFI: 235,014 / 46,995 / 24,862 across bins 1/2/3.
- Median σ comparison (`pooled`, `pos`, `neg`):
  - BMI: 0.00720 / 0.00610 / 0.00560
  - FIT: 0.00545 / 0.00512 / 0.00341
  - OFI: 0.04330 / 0.03210 / 0.03801
- Coverage validation: both `fm_stdev_pos_v3.RDS` and `fm_stdev_neg_v3.RDS` contain all 27 required table cells for BMI/FIT/OFI × spec_idx 1–3 × raw bins/differences, with no missing combinations in the validated subset.

### v3 vs v2 Notes
- Negative OFI spec3 changes only modestly from v2 to v3: bin1 +0.09, bin2 +0.07, bin3 +0.08.
- Negative FIT spec3 changes materially in bin1: v2 bin1 = 2.74 vs v3 bin1 = 0.80; bins 2/3 are similar or slightly larger. This reinforces that FIT sell-side results are sensitive to bin definition and should be treated as robustness evidence rather than the core pattern.

### Output
- `code/g_asym/archive/fm_stdev_pos_v3.RDS`
- `code/g_asym/archive/fm_stdev_neg_v3.RDS`
- `code/g_asym/archive/fm_nonlinear_pos_v3.RDS`
- `code/g_asym/archive/fm_nonlinear_neg_v3.RDS`
- `code/g_asym/archive/reg_asym_stdev_v3.tex` — generated; two-panel LaTeX table with raw multipliers, controls/stats rows, and differences block

### Direct-Mode Self-Review
- Describe-before-transform requirement satisfied: v3 prints sign-filter row counts, by-type distribution summaries, recomputed bin counts, and by-type sign-specific σ summaries before running regressions.
- Analyze/validate requirements satisfied: the sign filter logs retained row counts; bin construction is a single explicit transformation; output coverage checks found 0 missing required cells and 0 NaN coef/se values for both panels.
- Reproducibility verified by fresh `Rscript 2a_regression_fm_asym_v3.R` run before archival.

---

## Layout Cleanup

**Status:** IMPLEMENTED

### Current Layout
- Active asymmetry code and handoff docs now live in `code/g_asym/`; prior v1/v3 code and historical generated outputs live in `code/g_asym/archive/`.
- Upstream inputs remain in their existing locations: v1 RDS inputs in `20250117_quarterly/`; v2/v3 regression inputs and controls in `code/R/tmp/...`.
- `2a_regression_fm_asym_v2.R` sources utilities from `code/R/utilities/`, reads inputs from `code/R/tmp/...`, and writes generated v2 outputs into `code/g_asym/`.
- The active v2 note in `code/g_asym/asym_stdev_v2_note.tex` inputs the archived v2 table at `archive/reg_asym_stdev_v2.tex`.
- `code/README.md` now documents the `g_asym` folder.

### Verification
- `Rscript 2a_regression_fm_asym_v2.R` run from `code/g_asym/` successfully; wrote local v2 RDS outputs and `reg_asym_stdev_v2.tex`.
- `Rscript 2a_regression_fm_asym_v3.R` run from `code/g_asym/` successfully; wrote local v3 RDS outputs and `reg_asym_stdev_v3.tex`.
- `pdflatex -interaction=nonstopmode asym_stdev_v2_note.tex` run from `code/g_asym/` successfully; wrote `asym_stdev_v2_note.pdf`.
- Final RDS coverage check found 0 missing required cells and 0 NaN coef/se values for v2/v3 × pos/neg.
- `python3 table_asym_stdev.py` was not rerun because the current Python environment lacks `pyreadr`; the already-generated v1 table is present at `code/g_asym/archive/reg_asym_stdev.tex`.

### Current Cleanup State
- Top-level `code/g_asym/` intentionally keeps the active v2 script, active handoff docs, and the active v2 note.
- Historical outputs and non-v2 scripts are under `code/g_asym/archive/`.
- The active v2 script now prints sign-filter row counts, by-type `ofi` summaries, and by-bin count/mean-absolute-demand diagnostics before estimation.

---

## Task 4: Pooled Positive/Negative Bin Regression (v3)

**Status:** IMPLEMENTED and run

### Specification
- Script: `code/g_asym/2a_regression_fm_asym_v3.R`.
- Input: `code/R/tmp/raw_data/reg_inputs/reg_table_static.RDS`.
- Constructs six bin regressors for each demand type: `ofi_bin{k}_pos = ofi_bin{k} * 1{ofi > 0}` and `ofi_bin{k}_neg = ofi_bin{k} * 1{ofi < 0}`.
- Runs one pooled monthly cross-sectional regression per demand type/specification, with all six sign-bin variables included together.
- Keeps the v2 architecture: BMI/FIT/OFI all run; spec_idx 1-3 are no controls / return-predictor controls / return-predictor + liquidity controls; progressive interaction specs are also run.
- Uses the same v2 controls. BMI additionally receives the BMI-specific controls loaded from `controls_for_BMI.RDS`.

### Outputs
- `code/g_asym/fm_stdev_posneg_v3.RDS`
- `code/g_asym/fm_stdev_posneg_v3.dta`
- `code/g_asym/reg_asym_stdev_posneg_v3.tex`
- `code/g_asym/reg_asym_stdev_posneg_v3_balanced.tex`
- `code/g_asym/asym_stdev_v3_note.tex`
- `code/g_asym/asym_stdev_v3_note.pdf`
- `code/g_asym/make_asym_stdev_v3_note_table.R`
- `code/g_asym/make_asym_stdev_v3_balanced_note_table.R`

### Most-Controlled Results (spec_idx = 3)
- FIT positive bins: 6.95 / 5.39 / 3.65; differences 2-1 = -1.56 (t = -2.77), 3-2 = -1.74 (t = -3.06), 3-1 = -3.31 (t = -3.93).
- FIT negative bins: 3.58 / 2.17 / 2.25; differences 2-1 = -1.32 (t = -1.66), 3-2 = 0.27 (t = 0.38), 3-1 = -0.88 (t = -0.88).
- OFI positive bins: 3.66 / 2.49 / 1.69; all pairwise differences negative and large in t-stat magnitude.
- OFI negative bins: 2.67 / 1.71 / 0.60; all pairwise differences negative and large in t-stat magnitude.
- BMI remains noisy, with only 21 months in the BMI sample.

### Verification
- `Rscript 2a_regression_fm_asym_v3.R` ran successfully from `code/g_asym/`.
- Sign/bin diagnostics confirm mean absolute demand rises from bin 1 to bin 3 within each type and sign.
- Omitted monthly coefficients are retained as missing, not replaced with zero. This is intentional and differs from the Stata `statsby _b[var]` output, which effectively stores zero for omitted coefficients.
- FIT/spec3 omitted coefficients:
  - `ofi_bin2_neg`: 200106, 200203.
  - `ofi_bin3_neg`: 200106, 200203, 200403, 201303.
- These omissions occur because the corresponding negative FIT bin has zero observations in those months: 200106 has 0 bin2 / 0 bin3 observations; 200203 has 0 / 0; 200403 has 30 / 0; 201303 has 6 / 0.
- With missing omitted coefficients, FIT/spec3 negative-side summaries are: `ofi_bin2_neg` = 2.17 (SE 0.62, t 3.52, 117 months), `ofi_bin3_neg` = 2.25 (SE 0.56, t 3.98, 115 months), `d_neg_2_1` = -1.32 (t -1.66), `d_neg_3_2` = 0.27 (t 0.38), `d_neg_3_1` = -0.88 (t -0.88).
- The apparent mismatch between FIT/spec3 negative `ofi_bin2_neg` = 2.17, `ofi_bin3_neg` = 2.25, and `ofi_bin3_neg - ofi_bin2_neg` = 0.27 is due to different monthly samples. The coefficient rows average each coefficient over its own nonmissing months. The difference row averages the month-by-month paired difference over months where both coefficients are estimated. In the common-month sample, mean `ofi_bin2_neg` = 1.98 and mean `ofi_bin3_neg` = 2.25, so the paired difference is 0.27.
- If omitted coefficients are instead set to zero to mimic Stata's displayed `statsby` behavior, those values become `ofi_bin2_neg` = 2.14, `ofi_bin3_neg` = 2.17, `d_neg_2_1` = -1.45, `d_neg_3_2` = 0.04, and `d_neg_3_1` = -1.41. This zero-fill version is not the maintained R behavior.
- The v3 PDF note was later extended to include the balanced-month v3 variant table. It was recompiled with `pdflatex -interaction=nonstopmode asym_stdev_v3_note.tex`; the final run completed cleanly and wrote a 4-page PDF.

---

## Task 5: Rebuilt Positive/Negative Bin Regression (v4)

**Status:** IMPLEMENTED and run

### Specification
- Script: `code/g_asym/2a_regression_fm_asym_v4.R`.
- Input: `code/R/tmp/raw_data/reg_inputs/reg_table_static.RDS`.
- Keeps the v2 architecture and controls: BMI/FIT/OFI are run separately; spec_idx 1-3 are no controls / return-predictor controls / return-predictor + liquidity controls; progressive interaction specs are also run; BMI receives the BMI-specific controls from `controls_for_BMI.RDS`.
- Runs one pooled monthly cross-sectional regression with six sign-bin regressors.
- Rebuilds bins from raw `ofi` rather than using prebuilt `ofi_bin1/2/3`:
  - `ofi_sd_pos = sd(ofi | ofi > 0)` within `type, yyyymm`.
  - `ofi_sd_neg = sd(ofi | ofi < 0)` within `type, yyyymm`.
  - Positive bins use `ofi > 0`, `[0, sd_pos)`, `[sd_pos, 2 sd_pos)`, and `[2 sd_pos, inf)`.
  - Negative bins use `ofi < 0`, `(-sd_neg, 0)`, `[-2 sd_neg, -sd_neg]`, and `(-inf, -2 sd_neg]`.
- The `type, yyyymm` grouping is required because the all-type R script does not first `keep if type == ...`; it is equivalent to Stata's `bysort yyyymm` after filtering to one type.

### Outputs
- `code/g_asym/fm_stdev_posneg_rebuilt_v4.RDS`
- `code/g_asym/fm_stdev_posneg_rebuilt_v4.dta`
- `code/g_asym/reg_asym_stdev_posneg_rebuilt_v4.tex`
- `code/g_asym/summary_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex`
- `code/g_asym/reg_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex`
- `code/g_asym/asym_stdev_v4_note.tex`
- `code/g_asym/asym_stdev_v4_note.pdf`
- `code/g_asym/make_asym_stdev_v4_note_table.R`
- `code/g_asym/make_asym_stdev_v4_dynamic_summary_table.R`
- `code/g_asym/make_asym_stdev_v4_dynamic_note_table.R`

### Validation
- No missing sign-specific SD months:
  - BMI: 21 months, 0 missing positive-SD months, 0 missing negative-SD months.
  - FIT: 119 months, 0 missing positive-SD months, 0 missing negative-SD months.
  - OFI: 119 months, 0 missing positive-SD months, 0 missing negative-SD months.
- Bin counts match the rebuilt sign-specific v3 counts: BMI positive 3,765 / 670 / 584 and negative 3,671 / 477 / 533; FIT positive 223,728 / 52,235 / 31,910 and negative 155,449 / 49,096 / 30,160; OFI positive 163,014 / 37,776 / 20,083 and negative 235,014 / 46,995 / 24,862.
- Omitted monthly coefficients remain missing, not zero-filled.

### Most-Controlled Results (spec_idx = 3)
- BMI positive bins: 1.46 / 0.62 / 0.69; differences are noisy.
- BMI negative bins: 3.14 / 1.28 / 1.08; bin3-bin1 = -2.06 (t = -1.44).
- FIT positive bins: 7.32 / 6.09 / 4.14; bin3-bin1 = -3.18 (t = -3.13).
- FIT negative bins: 3.60 / 3.32 / 3.08; all pairwise differences are small and insignificant.
- OFI positive bins: 4.34 / 3.06 / 1.92; all pairwise differences are negative with large t-stat magnitudes.
- OFI negative bins: 2.93 / 1.85 / 0.70; all pairwise differences are negative with large t-stat magnitudes.
- The v4 PDF note now uses the v2 table layout: coefficient block, controls/Obs/R2 block, and pairwise-difference block inside each positive/negative panel. It also includes the dynamic v4 4-lag FIT/OFI summary-statistics table, the dynamic v4 4-lag FIT/OFI regression table, the most-controlled static v4 plots for BMI/FIT/OFI, and the most-controlled dynamic v4 4-lag plots for FIT/OFI. It was compiled with `pdflatex -interaction=nonstopmode asym_stdev_v4_note.tex`; the final run completed cleanly and wrote a 6-page PDF.

---

## Task 6: Balanced-Month v3 Pooled Positive/Negative Bin Regression

**Status:** IMPLEMENTED and run

### Specification
- Script: `code/g_asym/2a_regression_fm_asym_v3_balanced.R`.
- Input: `code/R/tmp/raw_data/reg_inputs/reg_table_static.RDS`.
- Same six-variable pooled v3 specification as `2a_regression_fm_asym_v3.R`.
- Additional sample restriction: within each demand type, drop `yyyymm` months where any of the six sign-bin regressors has zero nonzero observations. This makes the coefficient rows and paired-difference rows use the same month set for the sign-bin variables.
- Omitted coefficients are still kept missing if they occur for another reason; the filter is designed to eliminate the empty-sign-bin cases that drove the FIT negative mismatch in the unbalanced v3 table.

### Dropped Dates
- BMI: 21 total months, 21 kept, 0 dropped.
- FIT: 119 total months, 115 kept, 4 dropped.
- OFI: 119 total months, 119 kept, 0 dropped.
- Dropped FIT months:
  - 200106: `ofi_bin2_neg` = 0 obs, `ofi_bin3_neg` = 0 obs.
  - 200203: `ofi_bin2_neg` = 0 obs, `ofi_bin3_neg` = 0 obs.
  - 200403: `ofi_bin2_neg` = 30 obs, `ofi_bin3_neg` = 0 obs.
  - 201303: `ofi_bin2_neg` = 6 obs, `ofi_bin3_neg` = 0 obs.
- Rows after the filter: BMI 9,908; FIT 523,528; OFI 527,744.

### Outputs
- `code/g_asym/fm_stdev_posneg_v3_balanced.RDS`
- `code/g_asym/fm_stdev_posneg_v3_balanced.dta`
- `code/g_asym/fm_stdev_posneg_v3_balanced_dropped_type_months.RDS`
- `code/g_asym/fm_stdev_posneg_v3_balanced_dropped_type_months.dta`
- `code/g_asym/reg_asym_stdev_posneg_v3_balanced.tex`
- Included in `code/g_asym/asym_stdev_v3_note.pdf`.

### Most-Controlled Results (spec_idx = 3)
- FIT positive bins: 7.15 / 5.52 / 3.74; differences 2-1 = -1.62 (t = -2.80), 3-2 = -1.78 (t = -3.03), 3-1 = -3.41 (t = -3.93).
- FIT negative bins: 3.12 / 1.98 / 2.25; differences 2-1 = -1.15 (t = -1.47), 3-2 = 0.27 (t = 0.38), 3-1 = -0.88 (t = -0.88).
- OFI and BMI are unchanged relative to unbalanced v3 because no months are dropped for those demand types.

---

## Task 7: Presentation Plots for Pooled Asymmetry Results

**Status:** IMPLEMENTED and run

### Specification
- Script: `code/g_asym/pres_plots_static_asym_pooled.py`.
- Run command: `python3.11 pres_plots_static_asym_pooled.py`.
- The script mirrors the archived plotting format in `archive/pres_plots_static_asym_new.py`: one PNG per demand type, three x-axis comparisons (`Medium - Small`, `Big - Medium`, `Big - Small`), black positive-shock differences, red negative-shock differences, dashed 95% CIs, point labels, and bottom labels showing the relevant bin-coefficient subtraction.
- Uses the most-controlled specification (`var_added == "controls_char+controls_liq"`).
- Uses `python3.11` explicitly because the system `python3` is Python 3.9 with incompatible NumPy/pandas/matplotlib packages.

### Outputs
- v3 original pooled bins:
  - `code/g_asym/plots/v3/asym_v3_OFI.png`
  - `code/g_asym/plots/v3/asym_v3_FIT.png`
  - `code/g_asym/plots/v3/asym_v3_BMI.png`
- v3 balanced-month pooled bins:
  - `code/g_asym/plots/v3_balanced/asym_v3_balanced_OFI.png`
  - `code/g_asym/plots/v3_balanced/asym_v3_balanced_FIT.png`
  - `code/g_asym/plots/v3_balanced/asym_v3_balanced_BMI.png`
- v4 rebuilt bins:
  - `code/g_asym/plots/v4/asym_v4_OFI.png`
  - `code/g_asym/plots/v4/asym_v4_FIT.png`
  - `code/g_asym/plots/v4/asym_v4_BMI.png`

### Validation
- `python3.11` imports `pandas`, `matplotlib`, `numpy`, and `pyreadr` successfully.
- The script produced all 9 expected PNGs.
- Spot-checked `plots/v3_balanced/asym_v3_balanced_FIT.png`; it matches the archived visual format.

---

## Task 8: Dynamic Pooled Positive/Negative Bin Regression (v3)

**Status:** IMPLEMENTED and run

### Specification
- Script: `code/g_asym/2a_regression_fm_asym_v3_dynamic.R`.
- Input: `code/R/tmp/raw_data/reg_inputs/reg_table_dynamic.RDS`.
- Dynamic reference script: `code/R/c_dynamic_results/2_regression_dynamic_fm.R` (the prompt path `code/R/c_dynamic_results.R` is not present on disk).
- Mirrors the dynamic setup:
  - Melts `cumofi_1` through `cumofi_4` into `cumofi_lag`.
  - Creates type-horizon cells like `FIT_1lag`, `OFI_4lag`, and `OFI_pre_whitened_4lag`.
  - Bins are based on `abs(cumofi_lag)` within `yyyymm, type`.
  - The regressor value inside each bin is contemporaneous `ofi`, matching the existing dynamic script.
- Positive/negative sign is based on lagged cumulative shock, `cumofi_lag`, not contemporaneous `ofi`.
- Runs one pooled monthly cross-sectional regression with six variables:
  `ofi_bin1_pos`, `ofi_bin1_neg`, `ofi_bin2_pos`, `ofi_bin2_neg`, `ofi_bin3_pos`, `ofi_bin3_neg`.
- Uses the dynamic script's controls and progressive `ofi × control` interaction architecture.
- Produces both unbalanced and balanced variants:
  - Unbalanced keeps all type-horizon-months and leaves omitted monthly coefficients missing.
  - Balanced drops type-horizon-months where any of the six sign-bin regressors has zero nonzero observations.

### Input Diagnostics
- Dynamic rows after melting and control merge: 4,974,241.
- Type-horizon coverage:
  - FIT: 118 / 117 / 116 / 115 months for horizons 1--4.
  - OFI: 118 / 117 / 116 / 115 months for horizons 1--4.
  - OFI_pre_whitened: 116 / 115 / 114 / 113 months for horizons 1--4.
- Bin diagnostics confirm mean absolute `cumofi_lag` rises monotonically from bin 1 to bin 3 within demand type, horizon, and sign.

### Balanced Dropped Dates
- FIT_1lag: 4 months dropped — 200109, 200206, 200406, 201306.
- FIT_2lag: 7 months dropped — 200109, 200112, 200209, 200312, 200403, 200406, 201309.
- FIT_3lag: 8 months dropped — 200109, 200112, 200203, 200206, 200209, 200403, 200406, 201312.
- FIT_4lag: 7 months dropped — 200112, 200203, 200206, 200209, 200403, 200406, 200409.
- OFI_1lag and OFI_2lag: 0 months dropped.
- OFI_3lag: 3 months dropped — 201906, 201909, 202212.
- OFI_4lag: 5 months dropped — 200309, 201812, 201906, 201909, 202212.
- OFI_pre_whitened horizons 1--4: 0 months dropped.

### Outputs
- `code/g_asym/fm_dynamic_stdev_posneg_v3_unbalanced.RDS`
- `code/g_asym/fm_dynamic_stdev_posneg_v3_unbalanced.dta`
- `code/g_asym/fm_dynamic_stdev_posneg_v3_balanced.RDS`
- `code/g_asym/fm_dynamic_stdev_posneg_v3_balanced.dta`
- `code/g_asym/fm_dynamic_stdev_posneg_v3_balanced_dropped_type_months.RDS`
- `code/g_asym/fm_dynamic_stdev_posneg_v3_balanced_dropped_type_months.dta`

### 4-Lag Most-Controlled Results (spec_idx = 3)
- Unbalanced FIT positive bins: 3.85 / 3.52 / 1.55; differences 2-1 = -0.33 (t = -1.03), 3-2 = -1.98 (t = -3.87), 3-1 = -2.30 (t = -4.44).
- Unbalanced FIT negative bins: 3.65 / 1.76 / -2.33; differences 2-1 = -1.90 (t = -2.21), 3-2 = -4.84 (t = -1.39), 3-1 = -6.00 (t = -1.73).
- Unbalanced OFI positive bins: 1.51 / 0.63 / 0.60; differences 2-1 = -0.88 (t = -4.91), 3-1 = -0.94 (t = -5.18).
- Unbalanced OFI negative bins: 1.92 / 1.64 / 1.15; differences 2-1 = -0.29 (t = -3.40), 3-2 = -0.49 (t = -3.54), 3-1 = -0.77 (t = -5.92).
- Balanced FIT positive bins: 3.83 / 3.50 / 1.51; differences 2-1 = -0.32 (t = -0.97), 3-2 = -1.99 (t = -3.66), 3-1 = -2.31 (t = -4.20).
- Balanced FIT negative bins: 3.68 / 2.52 / -2.33; differences 2-1 = -1.16 (t = -2.30), 3-2 = -4.84 (t = -1.39), 3-1 = -6.00 (t = -1.73).
- Balanced OFI positive bins: 1.54 / 0.70 / 0.61; differences 2-1 = -0.85 (t = -4.66), 3-1 = -0.94 (t = -5.13).
- Balanced OFI negative bins: 1.93 / 1.64 / 1.16; differences 2-1 = -0.30 (t = -3.44), 3-2 = -0.48 (t = -3.38), 3-1 = -0.78 (t = -5.74).

---

## Task 9: Dynamic Rebuilt Positive/Negative Bin Regression (v4)

**Status:** IMPLEMENTED and run

### Specification
- Script: `code/g_asym/2a_regression_fm_asym_v4_dynamic.R`.
- Input: `code/R/tmp/raw_data/reg_inputs/reg_table_dynamic.RDS`.
- Same dynamic architecture as Task 8, but bins are rebuilt separately for positive and negative lagged cumulative shocks:
  - `cumofi_sd_pos = sd(cumofi_lag | cumofi_lag > 0)` within `yyyymm, type`.
  - `cumofi_sd_neg = sd(cumofi_lag | cumofi_lag < 0)` within `yyyymm, type`.
  - Positive and negative bin membership are based on `cumofi_lag`.
  - The regressor value inside each bin remains contemporaneous `ofi`.
- Produces unbalanced and balanced variants. The balanced variant drops type-horizon-months where any of the six sign-bin regressors has zero nonzero observations.

### Validation
- No missing positive- or negative-side SD months for any demand type or horizon.
- Rebuilt bin diagnostics show mean absolute `cumofi_lag` rises from bin 1 to bin 3 within each demand type, horizon, and sign.
- Balanced dynamic v4 drops 0 months in every type-horizon cell, so balanced and unbalanced outputs use the same samples.

### Outputs
- `code/g_asym/fm_dynamic_stdev_posneg_rebuilt_v4_unbalanced.RDS`
- `code/g_asym/fm_dynamic_stdev_posneg_rebuilt_v4_unbalanced.dta`
- `code/g_asym/fm_dynamic_stdev_posneg_rebuilt_v4_balanced.RDS`
- `code/g_asym/fm_dynamic_stdev_posneg_rebuilt_v4_balanced.dta`
- `code/g_asym/fm_dynamic_stdev_posneg_rebuilt_v4_balanced_dropped_type_months.RDS`
- `code/g_asym/fm_dynamic_stdev_posneg_rebuilt_v4_balanced_dropped_type_months.dta`
- `code/g_asym/summary_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex`
- `code/g_asym/reg_dynamic_asym_stdev_posneg_rebuilt_v4_4lag.tex`
- `code/g_asym/make_asym_stdev_v4_dynamic_summary_table.R`
- `code/g_asym/make_asym_stdev_v4_dynamic_note_table.R`
- `code/g_asym/asym_stdev_v4_note.tex`
- `code/g_asym/asym_stdev_v4_note.pdf`

### 4-Lag Most-Controlled Results (spec_idx = 3)
- FIT positive bins: 3.97 / 3.33 / 2.29; differences 2-1 = -0.64 (t = -1.67), 3-2 = -1.04 (t = -2.46), 3-1 = -1.68 (t = -4.33).
- FIT negative bins: 3.97 / 3.42 / 1.54; differences 2-1 = -0.55 (t = -1.14), 3-2 = -1.88 (t = -3.83), 3-1 = -2.43 (t = -4.57).
- OFI positive bins: 1.63 / 1.09 / 0.57; differences 2-1 = -0.54 (t = -6.05), 3-2 = -0.52 (t = -4.99), 3-1 = -1.06 (t = -9.41).
- OFI negative bins: 1.95 / 1.79 / 1.35; differences 2-1 = -0.16 (t = -2.13), 3-2 = -0.45 (t = -5.47), 3-1 = -0.60 (t = -7.29).
- The v4 PDF note reports these 4-lag dynamic results in the same v2-style layout as the static v4 table. It also includes the 4-lag dynamic bin summary table. The summary table reports bin counts and mean/SD for both the lagged cumulative demand shock used for binning and the contemporaneous demand shock used as the regression value.
- The dynamic table uses the unbalanced v4 output; balanced is numerically identical because zero type-horizon-months are dropped by the balanced filter. The note also includes the most-controlled dynamic v4 4-lag FIT/OFI plots.

---

## Task 10: Presentation Plots for Dynamic Pooled Asymmetry Results

**Status:** IMPLEMENTED and run

### Specification
- Script: `code/g_asym/pres_plots_dynamic_asym_pooled.py`.
- Run command: `python3.11 pres_plots_dynamic_asym_pooled.py`.
- Mirrors the archived dynamic/static plot format: one PNG per demand type and horizon, three x-axis comparisons (`Medium - Small`, `Big - Medium`, `Big - Small`), black positive-lagged-shock differences, red negative-lagged-shock differences, dashed 95% CIs, point labels, and bottom labels showing bin-coefficient subtraction.
- Uses the most-controlled specification (`var_added == "controls_char+controls_liq"`).
- Uses `python3.11` explicitly.

### Outputs
- 48 PNGs written under `code/g_asym/plots_dynamic/`:
  - `plots_dynamic/v3_unbalanced/`
  - `plots_dynamic/v3_balanced/`
  - `plots_dynamic/v4_unbalanced/`
  - `plots_dynamic/v4_balanced/`
- Each variant folder contains FIT, OFI, and OFI_pre_whitened plots for horizons 1--4.

### Validation
- The script completed successfully and printed all 48 output paths.
- `find plots_dynamic -name '*.png' | wc -l` returns 48.
- Spot-checked `plots_dynamic/v4_unbalanced/dynamic_asym_v4_FIT_4lag.png`; it matches the archived plot format.
