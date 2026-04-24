####################################################################################
# IMPORTS
####################################################################################

import pandas as pd
import pickle
import glob
import os
import matplotlib.pyplot as plt
import numpy as np
import pyreadr

from scipy import stats
import sys
from copy import deepcopy
from pprint import pprint
import re


####################################################################################
# LOAD GMM ESTIMATION RESULTS FROM ALL QUARTERLY WINDOWS
####################################################################################
# Each .p file contains the GMM estimation output from optimal_giv_FINAL.py for one
# rolling 4-quarter window.  The pickle stores:
#   all_results[0] = best solution tuple (zeta, info_dict, ier, ...)
#   all_results[1] = (unused here)
#   all_results[2] = list of ALL starting-point solutions: [(zeta, info_dict, ier), ...]
#   all_results[3] = (unused here) 
#   all_results[4] = starting points array
#
# For each window we select the best root via the J-statistic (Appendix C.2.2).
####################################################################################

# outer_folder = 'HOLDINGS_ESTIMATION_RESULTS/numFactors_5/'
# outer_folder = 'HOLDINGS_ESTIMATION_RESULTS_OSC/numFactors_5/'
outer_folder = 'HOLDINGS_ESTIMATION_RESULTS_OSC2/numFactors_5/'



all_files = os.listdir(outer_folder)
p_files = [f for f in all_files if f.endswith('.p')]
list_of_filenames = sorted([os.path.join(outer_folder, f) for f in p_files])
characteristics = ['active_share']

results_list = []


for filename in list_of_filenames:
 
    year_quarter = filename[filename.rfind('/') + 1 : filename.rfind('.')]
    year = int(year_quarter[:4])

    all_results = pickle.load(open(filename, 'rb'))
    for s_index, s in enumerate(all_results[2]):
        sse = sum(s[1]['fvec'] ** 2)
        all_results[2][s_index][1]['sse'] = sse

        

    # Step 1: Among converged roots with zeta_{1,0} > 0, pick the one minimizing
    #         the over-identified J-statistic (eq. C.13).  The economic prior
    #         zeta_{1,0} > 0 means demand curves slope downward.
    all_converged_roots = sorted([ (s[1]['J_stats']['moment_conditions_fast_overIdentified'], s)  for s in all_results[2] if s[2] == 1 and s[0][0] > 0], key = lambda s : s[0])

    # Step 2: Fall back to all converged roots (regardless of sign)
    if len(all_converged_roots) <= 0:
        all_converged_roots =  sorted([ (s[1]['J_stats']['moment_conditions_fast_overIdentified'], s)  for s in all_results[2] if s[2] == 1], key = lambda s : s[0])

    # Step 3: Fall back to all roots (including non-converged)
    if len(all_converged_roots) <= 0:
        all_converged_roots =  sorted([ (s[1]['J_stats']['moment_conditions_fast'], s)  for s in all_results[2] ], key = lambda s : s[0])


    result_fsolve = all_converged_roots[0][1]
    
    results_list.append([ pd.to_datetime(filename[filename.rfind('/') + 1 : filename.find('.') ]) ] + list(result_fsolve[0]))

    results_list[-1].append(result_fsolve[1]['nobs'])


####################################################################################
# ASSEMBLE RESULTS INTO A TIME-SERIES DATAFRAME
####################################################################################
# Each row = one rolling window (indexed by ending quarter).
# Columns = [zeta_{1,0}, zeta_{1,AS}, zeta_{2,0}, zeta_{2,AS}, nobs]
# where "NL" (nonlinear) corresponds to zeta_2 parameters.
####################################################################################

results_df = pd.DataFrame(results_list)

char_names = ['intercept'] + list(characteristics)
param_cols = [f"Param_{char}" for char in char_names ] + [f"Param_NL_{char}" for char in char_names]

results_df.columns = ['Date'] + param_cols + ['nobs']
    
results_df = results_df.set_index('Date')
 

####################################################################################
# FAMA-MACBETH STANDARD ERRORS WITH NEWEY-WEST CORRECTION
####################################################################################
# Table 6 reports time-series average estimates with Fama-MacBeth standard errors
# adjusted for autocorrelation using a Newey-West kernel (Section 5.1.3).
# Overlapping 4-quarter rolling windows induce serial correlation, so
# we use 8 lags (two full estimation windows) for the Newey-West correction.
####################################################################################

def get_autocorr_adj_FM_ses(results_df):
    """
    Compute Newey-West autocorrelation-adjusted Fama-MacBeth standard errors.

    Inputs:
        results_df : pd.DataFrame (T_windows x P), time series of parameter
                     estimates from each rolling window, indexed by Date.
                     Columns include 'Param_intercept', 'Param_active_share',
                     'Param_NL_intercept', 'Param_NL_active_share', 'nobs'.

    Output:
        pd.DataFrame (P x 3) with columns ['Coeff', 'SE', 'T'] where
            Coeff = time-series mean of each parameter
            SE    = Newey-West standard error of that mean
            T     = t-statistic = Coeff / SE
    """
    max_lag = 8

    print('max_lag: ' + str(max_lag))

    col_list = [c for c in  results_df.columns if 'Param'in c]

    coeffs = []
    ses = []
    t_stats = []

    for col in col_list:

        covs = np.zeros(max_lag + 1)

        numeric_series = pd.to_numeric(results_df[col], errors='coerce').dropna()
        mean = np.nanmean(numeric_series.values)

        T = numeric_series.shape[0]
        weights = 1 - np.arange(max_lag + 1 )/( max_lag + 1 )

        for l in range(max_lag + 1):
                
            cov_l = (numeric_series - mean) \
                    * (numeric_series.shift(l) - mean) \
                    / T
                    
            cov_l = cov_l.sum()

            weight =  weights[l]

            covs[l] = cov_l * weight
                
        var = 1 / numeric_series.shape[0] * ( covs[0] + 2 * covs[1:].sum() )
 
        coeffs.append(mean)
        ses.append(np.sqrt(var))
        t_stats.append( mean / np.sqrt(var) )

    df = pd.concat([pd.DataFrame(coeffs), pd.DataFrame(ses), pd.DataFrame(t_stats)], axis = 1)
    df.index = col_list
    df.columns = ['Coeff', 'SE', 'T']

    return df
  


def create_latex_table(df, within_q_lag_cum_ret_std):
    """
    Generate a LaTeX table of elasticity estimates (Table 6 in the paper).

    The table displays the four estimated parameters:
        zeta_{1,0}              : average elasticity (intercept)
        zeta_{1,Active Share}   : elasticity interaction with active share
        zeta_{2,0}              : sensitivity of elasticity to |P_tilde| (intercept)
        zeta_{2,Active Share}   : sensitivity interaction with active share

    Inputs:
        df                       : pd.DataFrame (P x 3), output from
                                   get_autocorr_adj_FM_ses with cols ['Coeff','SE','T']
        within_q_lag_cum_ret_std : float, cross-sectional std dev of
                                   |sum_l Delta_p_{n,t-l}| (from Table 5)

    Output:
        str, LaTeX source for the table
    """
    
    def format_num(x, decimals=4):
        return f"{x:.{decimals}f}"
    
    def add_stars(coef, se):
        z_score = abs(float(coef)/float(se))
        p_value = 2 * (1 - stats.norm.cdf(z_score))
        stars = ''
        if p_value < 0.01:
            stars = '^{***}'
        elif p_value < 0.05:
            stars = '^{**}'
        elif p_value < 0.1:
            stars = '^{*}'
        return f"${format_num(float(coef))}{stars}$"
    
    latex = "\\begin{tabular}{lcccc}\n"
    latex += "\\toprule\n"
    
    latex += " & $\zeta_{1,0}$ & $\zeta_{1,\\text{Active Share}}$ & $\zeta_{2,0}$ & $\zeta_{2,\\text{Active Share}}$ \\\\\n"

    latex += "\\midrule\n"
    
    params = ['Param_intercept', 'Param_active_share', 'Param_NL_intercept', 'Param_NL_active_share']
    coeff_values = []
    se_values = []
    
    for param in params:
        if param in df.index:
            coeff = df.loc[param, 'Coeff']
            se = df.loc[param, 'SE']
            coeff_values.append(add_stars(coeff, se))
            se_values.append(f"({format_num(se)})")
        else:
            coeff_values.append("--")
            se_values.append("--")
    
    latex += "Coefficient & " + " & ".join(coeff_values) + "\\\\\n"
    
    latex += " & " + " & ".join(se_values) + "\\\\\n"
    
    latex += "\\bottomrule\n"
    latex += "\\end{tabular}"
    
    return latex

####################################################################################
# COMPUTE TABLE 6: DEMAND ELASTICITY PARAMETER ESTIMATES
####################################################################################
# Rescaling:
#   - Active share coefficients are multiplied by 2 because the raw active share
#     has range [0, 0.5] but is reported on a [0, 1] scale in the paper.
#   - zeta_2 coefficients are multiplied by 100 because cumulative returns are
#     stored in percentage points internally but the paper reports them in
#     percentage-point units (0-100 scale).
####################################################################################

df = get_autocorr_adj_FM_ses(results_df)
df.loc['Param_active_share', 'Coeff'] *= 2
df.loc['Param_active_share', 'SE'] *= 2
df.loc['Param_NL_active_share', 'Coeff'] *= 2 * 100
df.loc['Param_NL_active_share', 'SE'] *= 2 * 100

df.loc['Param_NL_intercept', 'Coeff'] *= 100
df.loc['Param_NL_intercept', 'SE'] *= 100

within_q_lag_cum_ret_std = 32.25	
 
latex_table = create_latex_table(df, within_q_lag_cum_ret_std)
print(latex_table)

with open('giv_results_table.tex', 'w') as f:
    f.write(latex_table)



 
####################################################################################
# FIGURE 2: VISUALIZE VARIATION IN PRICE ELASTICITIES
####################################################################################
# Figure 2 in the paper visualizes how total elasticity
#   zeta(a, P_tilde) = zeta_{1,0} + zeta_{1,AS} * (a - mean_a)
#                     + zeta_{2,0} * (P_tilde - mean_P_tilde)
#                     + zeta_{2,AS} * (a - mean_a) * (P_tilde - mean_P_tilde)
# varies with active share (a) and absolute lagged cumulative price change
# (P_tilde = |sum_{l=1}^L Delta_p_{n,t-l}|).
#
# Panel (a): heatmap of elasticity over (P_tilde, active_share) grid
# Panel (b): elasticity vs P_tilde for fixed active share levels (p25, mean, p75)
####################################################################################

df = get_autocorr_adj_FM_ses(results_df)

# Extract the four zeta coefficients (raw estimation scale)
zeta_10, zeta_1_active_share, zeta_20, zeta_2_active_share = df['Coeff'].values

# Rescale active share coefficients to [0,1] scale
zeta_1_active_share *= 2
zeta_2_active_share *= 2

# Sample means from Table 5 (used to demean active share and P_tilde)
MEAN_ACTIVE_SHARE = 0.88 / 2
MEAN_LAG_RET = 31.55
active_share_grid_1d = np.linspace(0.20, 1.54, 150) / 2

 
def zeta(active_share, lag_ret, zeta_10, zeta_1_active_share, zeta_20, zeta_2_active_share, MEAN_ACTIVE_SHARE, MEAN_LAG_RET):
    return zeta_10 + zeta_1_active_share*(active_share - MEAN_ACTIVE_SHARE) + zeta_20*(lag_ret - MEAN_LAG_RET) + zeta_2_active_share*((active_share - MEAN_ACTIVE_SHARE)*(lag_ret - MEAN_LAG_RET))


# Grid over 5th-95th percentile ranges (from Table 5): P_tilde in [0, 81]
lag_ret_grid_1d = np.linspace(0	, 81, 150)

active_share_grid, lag_ret_grid = np.meshgrid(active_share_grid_1d, lag_ret_grid_1d, indexing="ij")
Z = zeta(active_share_grid, lag_ret_grid, zeta_10, zeta_1_active_share, zeta_20, zeta_2_active_share, MEAN_ACTIVE_SHARE, MEAN_LAG_RET) 

# ---- Panel (a): Heatmap of elasticity over (P_tilde, active_share) ----
plt.figure()
m = plt.pcolormesh(lag_ret_grid, active_share_grid, Z, shading="auto", cmap = 'RdYlBu')

contours = plt.contour(lag_ret_grid, active_share_grid, Z, colors='black', alpha=0.6, linewidths=2)

labels = plt.clabel(contours, inline=True, fontsize=12, fmt='%.2f')
for label in labels:
    label.set_color('black')
 

plt.xlabel("Absolute Lag Annual Price Change (%)", fontsize = 20)
plt.ylabel("Active Share", fontsize = 20)
cbar = plt.colorbar(m, label="Elasticity")
cbar.ax.tick_params(labelsize=20)
cbar.set_label("Elasticity", fontsize=20)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.axhline(y=MEAN_ACTIVE_SHARE, color='black', linestyle='--', linewidth=1, label = 'Mean Values')
plt.axvline(x=MEAN_LAG_RET, color='black', linestyle='--', linewidth=1)
plt.legend(fontsize=16, loc='best')

plt.gcf().set_size_inches(14, 10)
plt.savefig('heatmap_activeShare_vs_lagRet_elasticity.png', bbox_inches='tight')
plt.clf()

# ---- Panel (b): Elasticity vs P_tilde for fixed active shares (p25, mean, p75) ----

lag_ret_grid_1d_sub_grid = lag_ret_grid_1d


line_colors = ['#d62728', '#000000', '#1f77b4']

# Active share values from Table 5: 25th percentile, mean, 75th percentile
active_share_levels = {'25$^{th}$ Percentile': 0.54/2, 'Mean': 0.82/2, '75$^{th}$ Percentile': 1.10/2}
plt.figure()
for i, (label, a0) in enumerate(active_share_levels.items()):
    y = zeta(a0, lag_ret_grid_1d_sub_grid, zeta_10, zeta_1_active_share, zeta_20, zeta_2_active_share, MEAN_ACTIVE_SHARE, MEAN_LAG_RET)
    plt.plot(lag_ret_grid_1d_sub_grid, y, color=line_colors[i], label=f"Active Share = {label} ({a0:0.3f})", linewidth = 2)
    print(label, a0, y)

plt.axvline(MEAN_LAG_RET, linestyle="--", linewidth=1.5, color='black', label='Absolute Lag Annual Price Change Mean')

plt.xlabel("Absolute Lag Annual Price Change (%)", fontsize = 20)
plt.ylabel("Elasticity", fontsize = 20)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.legend(fontsize=16, loc='best')

plt.gcf().set_size_inches(14, 10)
plt.savefig('lageRet_elasticity_cuts.png', bbox_inches='tight')
plt.clf()



