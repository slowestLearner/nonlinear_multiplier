# Asymmetry Table — Results

> Mirrors PLAN.md structure. Updated after each step with key findings.
> New agents: read PLAN.md for what to do, RESULTS.md for what was found.

**Last updated:** 2026-04-26 (Task 3 v3 run + validation)
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
- `code/g_asym/reg_asym_stdev.tex` — 2,472 chars, two-panel LaTeX tabular (differences only; v1 data)

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
- Script `code/g_asym/2a_regression_fm_asym_v3.R` run successfully.
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
- `code/g_asym/fm_stdev_pos_v3.RDS`
- `code/g_asym/fm_stdev_neg_v3.RDS`
- `code/g_asym/fm_nonlinear_pos_v3.RDS`
- `code/g_asym/fm_nonlinear_neg_v3.RDS`
- `code/g_asym/reg_asym_stdev_v3.tex` — generated; two-panel LaTeX table with raw multipliers, controls/stats rows, and differences block

### Direct-Mode Self-Review
- Describe-before-transform requirement satisfied: v3 prints sign-filter row counts, by-type distribution summaries, recomputed bin counts, and by-type sign-specific σ summaries before running regressions.
- Analyze/validate requirements satisfied: the sign filter logs retained row counts; bin construction is a single explicit transformation; output coverage checks found 0 missing required cells and 0 NaN coef/se values for both panels.
- Reproducibility verified by fresh `Rscript 2a_regression_fm_asym_v3.R` run from `code/g_asym`.

---

## Layout Cleanup

**Status:** IMPLEMENTED

### Current Layout
- Active asymmetry code, handoff docs, generated LaTeX tables, note files, plots, and v2/v3 RDS outputs now live in `code/g_asym/`.
- Upstream inputs remain in their existing locations: v1 RDS inputs in `20250117_quarterly/`; v2/v3 regression inputs and controls in `code/R/tmp/...`.
- `2a_regression_fm_asym_v2.R` and `2a_regression_fm_asym_v3.R` now source utilities from `code/R/utilities/`, read inputs from `code/R/tmp/...`, and write all generated outputs into `code/g_asym/`.
- `table_asym_stdev.py` now writes `reg_asym_stdev.tex` into `code/g_asym/`.
- `code/README.md` now documents the `g_asym` folder.

### Verification
- `Rscript 2a_regression_fm_asym_v2.R` run from `code/g_asym/` successfully; wrote local v2 RDS outputs and `reg_asym_stdev_v2.tex`.
- `Rscript 2a_regression_fm_asym_v3.R` run from `code/g_asym/` successfully; wrote local v3 RDS outputs and `reg_asym_stdev_v3.tex`.
- `pdflatex -interaction=nonstopmode asym_stdev_v2_note.tex` run from `code/g_asym/` successfully; wrote `asym_stdev_v2_note.pdf`.
- Final RDS coverage check found 0 missing required cells and 0 NaN coef/se values for v2/v3 × pos/neg.
- `python3 table_asym_stdev.py` was not rerun because the current Python environment lacks `pyreadr`; the already-generated v1 table is present at `code/g_asym/reg_asym_stdev.tex`.
