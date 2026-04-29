import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
from matplotlib.patches import FancyBboxPatch
import pyreadr
import pandas as pd

import os
os.environ["PATH"] = "/usr/local/texlive/2020/bin/x86_64-darwin:" + os.environ["PATH"]
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath,amssymb}"


mpl.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r"""
\usepackage[utf8]{inputenc}
\usepackage{amsmath,amssymb}
\usepackage{newunicodechar}
\newunicodechar{Δ}{$\Delta$} % map Unicode Delta to math \Delta
\newunicodechar{−}{-}        % map Unicode minus (U+2212) to ASCII hyphen
""",
})

# Tell Matplotlib to use LaTeX, but with DejaVu Sans for everything
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "text.latex.preamble": r"""
\usepackage[utf8]{inputenc}
\usepackage{sfmath}       % use sans-serif fonts in math mode
\usepackage{newunicodechar}
\newunicodechar{Δ}{$\Delta$}
\newunicodechar{−}{-}
""",
})

# ========== Read data ==========
pos = pyreadr.read_r('fm_stdev_pos.RDS')
neg = pyreadr.read_r('fm_stdev_neg.RDS')

pos_df = pos[None]
neg_df = neg[None]

# Filter for controls_char+controls_liq
pos_filt = pos_df[pos_df['var_added'] == 'controls_char+controls_liq']
neg_filt = neg_df[neg_df['var_added'] == 'controls_char+controls_liq']

# ========== Helper functions ==========

def base_plot(ax, y_max):
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-y_max*1.6, y_max*0.3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(bins, fontsize=FS_TICK)
    ax.set_ylabel(r"Multiplier Difference", fontsize=FS_YLAB)
    ax.axhline(0, color='0.7', linestyle='-', linewidth=1.5, zorder=1)
    for s in ['top','right']:
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='y', labelsize=FS_TICK)

def draw_point_with_cis(ax, i, coef_val, se_val, color, offset=0, bin_labels=None):
    ci95 = 1.96 * se_val
    # 95% CI (dashed)
    ax.vlines(i+offset, coef_val-ci95, coef_val+ci95,
              color=color, ls='--', lw=LW_95, alpha=ALPHA)
    ax.hlines([coef_val-ci95, coef_val+ci95], i+offset-cap_hw, i+offset+cap_hw,
              color=color, ls='--', lw=LW_95, alpha=ALPHA)
    # point
    ax.plot(i+offset, coef_val, 'o', color=color, markersize=MS)
    # label at point
    ax.text(i+offset+0.08, coef_val, f"{coef_val:.1f}",
            ha='left', va='center', fontsize=FS_LAB, color=color)
    # label at bottom of CI with bin difference
    if bin_labels is not None:
        ax.text(i+offset, coef_val-ci95-0.2, bin_labels[i],
                ha='center', va='top', fontsize=FS_BOTTOM, color=color)

def add_errorbar_legend(ax):
    """Add legend explaining positive/negative and CI bars."""
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], color='black', marker='o', linestyle='None',
               markersize=MS, label='Positive'),
        Line2D([0], [0], color='red', marker='o', linestyle='None',
               markersize=MS, label='Negative'),
        Line2D([0], [0], color=(0, 0, 0, ALPHA), lw=LW_95, ls='--', label='95\% CI')
    ]
    leg = ax.legend(handles=legend_elems,
                    fontsize=FS_LEG,
                    loc='lower left',
                    frameon=True,
                    facecolor='white',
                    edgecolor='0.6',
                    framealpha=1.0,
                    handlelength=2.8,
                    bbox_to_anchor=(0.02, 0.02))
    for text in leg.get_texts():
        text.set_color('0.1')


# ========== Data & display settings ==========

bins = ['Medium - Small', 'Big - Medium', 'Big - Small']

# Figure aesthetics
FS_YLAB = 34
FS_TITLE = 40
FS_TICK = 30
FS_LAB  = 26
FS_ANN  = 30
FS_LEG  = 24
FS_BOTTOM = 20
MS = 16
LW_95 = 1.8
ALPHA= 0.35
cap_hw = 0.03

names_dict = {
    'BMI': 'Benchmarking Intensity (BMI)',
    'FIT': 'Flow Induced Trading (FIT)',
    'OFI': 'Order Flow Imbalance (OFI)',
}

# Offset for parallel points
OFFSET = 0.20

# ========== Loop through types ==========
for type_name in ['OFI', 'FIT', 'BMI']:
    # Extract data
    pos_type = pos_filt[pos_filt['type'] == type_name]
    neg_type = neg_filt[neg_filt['type'] == type_name]

    # Get bin values for labels
    bin_vars = ['ofi_bin1', 'ofi_bin2', 'ofi_bin3']
    pos_bins = []
    neg_bins = []
    for var in bin_vars:
        pos_bins.append(pos_type[pos_type['var'] == var]['coef'].values[0])
        neg_bins.append(neg_type[neg_type['var'] == var]['coef'].values[0])

    # Get coefficients and SEs for the three differences
    vars_list = ['ofi_bin2 - ofi_bin1', 'ofi_bin3 - ofi_bin2', 'ofi_bin3 - ofi_bin1']

    pos_coefs = []
    pos_ses = []
    neg_coefs = []
    neg_ses = []

    for var in vars_list:
        pos_row = pos_type[pos_type['var'] == var]
        neg_row = neg_type[neg_type['var'] == var]

        pos_coefs.append(pos_row['coef'].values[0])
        pos_ses.append(pos_row['se'].values[0])
        neg_coefs.append(neg_row['coef'].values[0])
        neg_ses.append(neg_row['se'].values[0])

    pos_coefs = np.array(pos_coefs)
    pos_ses = np.array(pos_ses)
    neg_coefs = np.array(neg_coefs)
    neg_ses = np.array(neg_ses)

    # Create labels showing actual multiplier differences
    pos_labels = [
        f"{pos_bins[1]:.1f} - {pos_bins[0]:.1f}",  # bin2 - bin1
        f"{pos_bins[2]:.1f} - {pos_bins[1]:.1f}",  # bin3 - bin2
        f"{pos_bins[2]:.1f} - {pos_bins[0]:.1f}"   # bin3 - bin1
    ]
    neg_labels = [
        f"{neg_bins[1]:.1f} - {neg_bins[0]:.1f}",  # bin2 - bin1
        f"{neg_bins[2]:.1f} - {neg_bins[1]:.1f}",  # bin3 - bin2
        f"{neg_bins[2]:.1f} - {neg_bins[0]:.1f}"   # bin3 - bin1
    ]

    # Calculate y-axis range
    all_vals = np.concatenate([pos_coefs - 1.96*pos_ses, pos_coefs + 1.96*pos_ses,
                               neg_coefs - 1.96*neg_ses, neg_coefs + 1.96*neg_ses])
    y_max = np.max(np.abs(all_vals))

    # Create plot
    fig, ax = plt.subplots(figsize=(14, 10))
    base_plot(ax, y_max)

    # Draw points for all three comparisons
    for i in range(3):
        draw_point_with_cis(ax, i, pos_coefs[i], pos_ses[i], 'black', -OFFSET, pos_labels)
        draw_point_with_cis(ax, i, neg_coefs[i], neg_ses[i], 'red', OFFSET, neg_labels)

    add_errorbar_legend(ax)
    plt.subplots_adjust(bottom=0.15)
    plt.title(names_dict[type_name], fontsize=FS_TITLE, pad=50)
    plt.savefig(f"asym_{type_name}.png", dpi=300, bbox_inches='tight')
    plt.close()

print("Generated plots: asym_OFI.png, asym_FIT.png, asym_BMI.png")
