# Nonlinear multipliers

This repo accompanies paper:
[Chaudhry andLi (2024), "Endogenous Elasticities: Price Multipliers Are Smaller for Larger Demand Shocks"](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=5325371)

If you use data or code based on our work, please cite the paper: 

~~~
@article{chaudhry2025endogenous,
  title={Endogenous elasticities: Price multipliers are smaller for larger demand shocks},
  author={Chaudhry, Aditya and Li, Jiacui},
  year={2025}
}
~~~



----

## Data 

TODO: fill in details later

TODO later
1. OFI needs to be filled to 2023 (start from direct download from WRDS)
2. 


----

## R Code 



-`0_not_to_share`: basic data processing, keep to ourselves
	- `1_copy_upstream_data.R`: DOING

TODO for the final version
- remove pre-whitened OFI
- remove any reference to "nonlinear specification"

 
-`a_data_processing`
	-DONE. `1_get_daily_ofi_residuals.R`: Use rolling regressions to extract the surprise component of OFI
	-DONE. `2_put_together_main_data.R`: Combines quarterly returns with FIT, OFI, and BMI shocks; applies BMI bandwidth filters; saves core regression inputs and BMI controls.
	  - changed: only have data from 1993 to 2022
	-DONE. `3_sum_stats_cross_section.R`: Computes time-series averages of cross-sectional moments for demand/return variables and saves summary tables.
	-DONE. `4a_construct_liquidity_controls.R`: 
	-DONE. `4b_put_together_controls.R`: Merges liquidity and characteristics controls at quarterly freq for later regressions. Also create control metadata labels.
  

- `b_static_results`

	-DONE. `1_construct_regression_table_static.R`: Prepares static regression dataset by merging demand/returns with controls and generating nonlinear and shock-bin regressors.
	-DONE. `2a_regression_fm.R`: Runs static Fama-MacBeth regressions (nonlinear and bin-based) across demand types and progressive control specifications.
	-DONE. `2b_regression_panel.R`: Runs two-way-clustered panel versions of static regressions (mainly bin-based specification) as a robustness counterpart.
	-DONE. `3a_standardized_compute_lagged_vol.R`: Computes rolling 4/8/12-quarter demand volatility measures used to standardize FIT/OFI shocks.
	-DONE. `3b_standardized_regression.R`: Standardizes demand by lagged volatility and reruns static Fama-MacBeth stdev-bin regressions by volatility horizon.
	-DONE. `4_anticipation.R`: Tests anticipation by regressing cumulative future returns (multiple lags) on current demand with lag-adjusted controls.
	-DONE. `5_substitution.R`: Re-estimates static regressions allowing control-by-bin interactions to study substitution/heterogeneous effects.

  
- `c_dynamic_results`
	-DONE. `1_construct_regression_table_dynamic.R`: Converts static inputs into dynamic-regression inputs by adding lagged and cumulative past-demand terms.
	-DONE. `2_regression_dynamic_fm.R`: Runs dynamic Fama-MacBeth regressions linking returns to lagged cumulative demand with nonlinear and bin-based specifications.
	-DONE. `3_regression_dynamic_and_static_fm.R`: Specification that jointly uses contemporaneous and lagged demand bins with custom FM comparisons.


- `d_additional`
	-DONE `1_bmi_pass_thru.R`: Estimates whether BMI shock bins explain cross-sectional variation in institutional ownership pass-through (`dio`) with richer controls.
	-DONE `2a` to `2b`: taken from tests/27_alternative_fit_construtions/: clean FIT. 
	  - Difference: we NO LONGER take out fund-level characteristics and we ONLY take out PCs
  -DONE `2c_clean_fit_summary.R`: Summarizes alternative FIT/flow constructions via binned cross-sectional comparisons
	-DONE `2d_clean_fit_regression_fm.R`: Replaces FIT with alternative flow-based FIT constructions and reruns static Fama-MacBeth bin regressions.
	-DONE`3_reversal.R`: Tests return reversal patterns by regressing multi-horizon future returns on demand bins and interaction-enriched controls.
	-DONE `4a_fit_cleaning_prepare_data.R`: 
	-DONE `4b_fit_cleaning_get_results.R`: - somehow results are gone? No longer see clear differences in flow->trade coefs? Seriously? 

- `e_produce_tables`
	-DONE. `reg_bmi_pass_thru.Rmd`: Builds a two-panel LaTeX table for BMI pass-through regressions, showing bin coefficients for `ΔIO` plus coefficient-difference tests across specifications.
	-DONE. `reg_contemp_and_dynamic_fm_stdev.Rmd`: Produces a table that jointly estimates contemporaneous and dynamic (lagged cumulative-demand) FM price-impact effects for FIT and OFI.
	-DONE. `reg_dynamic_fm_stdev.Rmd`: Produces the main dynamic FM stdev table using lagged cumulative-demand bins, with controls and coefficient-comparison panel.
	-DONE. `reg_dynamic_fm_stdev_vary_lookback.Rmd`: Outputs a dynamic FM stdev table that holds specification fixed and varies lookback windows to compare bin effects by lag horizon.
	-DONE. `reg_static_fm_stdev.Rmd`: Produces the baseline static FM stdev table for bin-based demand effects, including controls, fit metrics, and coefficient-difference tests.
	-DONE. `reg_static_fm_stdev_anticipation.Rmd`: Builds an anticipation table that regresses returns on lagged demand bins (`h=1,2,4`) to test pre-return demand signals.
	-DONE. `reg_static_fm_stdev_fit_variations.Rmd`: Creates a FIT-robustness table comparing static FM bin coefficients as increasing numbers of FIT principal components are removed.
	-DONE. `reg_static_fm_stdev_reversal.Rmd`: Produces a reversal-horizon table linking current demand bins to future returns (`h=1,2,4`) across BMI/FIT/OFI.
	-DONE. `reg_static_fm_stdev_with_demand_interactions.Rmd`: static price-impact specifications with demand interactions.
	-DONE. `reg_static_fm_stdev_standardized.Rmd`: Generates a standardized-demand table comparing FIT/OFI bin effects across standardized lookback variants versus unstandardized baselines.
	-DONE. `reg_static_fm_stdev_substitution.Rmd`: Creates a substitution-specification table that swaps in predictor/liquidity-bin interaction structures and reports corresponding coefficient contrasts.
	-DONE.  `reg_static_panel_stdev.Rmd`: Generates a static panel-regression counterpart to the FM table for stdev-binned demand effects and coefficient-comparison diagnostics.
	-DONE. `summary_statistics.Rmd`: Produces the summary-statistics LaTeX table (obs, moments, and percentiles) for returns, demand proxies, and market cap.

## `f_produce_figures`

- DONE `2_progressively_add_controls.R`: Plots how key multiplier spread estimates evolve as interaction controls are incrementally added
- DONE `3_flow_cleaning.R`: Produces diagnostic figures for flow cleaning choices (bin behavior, trade-flow response, winsorization, and FIT before/after cleaning).
- DONE `4_flow_PCA.R`: Plots how removing principal components from flows changes the FIT mapping relative to original FIT.

  

## `utilities`

- `regressions.R`: Defines reusable estimation helpers for Fama-MacBeth and panel regressions, including coefficient-difference and Newey-West logic.
- `runmefirst.R`: Clears workspace, loads shared libraries/options, and defines a default parallel-core count used by downstream scripts.

  

## `renv`

- `activate.R`: Standard `renv` project autoloader script that bootstraps and activates the project-specific package library.
