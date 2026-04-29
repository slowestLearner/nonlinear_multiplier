import matplotlib.pyplot as plt
import numpy as np
from matplotlib import transforms
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
from matplotlib.patches import FancyBboxPatch

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
    # Optional: consistent math fonts
    # "font.family": "serif",
    # "font.serif": ["Times New Roman"],
})

# (2) Tell Matplotlib to use LaTeX, but with DejaVu Sans for everything
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],  # <- this is Matplotlib’s normal default
    "text.latex.preamble": r"""
\usepackage[utf8]{inputenc}
\usepackage{sfmath}       % use sans-serif fonts in math mode
\usepackage{newunicodechar}
\newunicodechar{Δ}{$\Delta$}
\newunicodechar{−}{-}
""",
})

# ========== Helper functions ==========

def base_plot(ax):
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-0.2, max(coef + ci95) + 0.8)
    ax.set_xticks(range(3))
    ax.set_xticklabels(bins, fontsize=FS_TICK)
    ax.set_ylabel(r"Price Multiplier", fontsize=FS_YLAB)
    for s in ['top','right']:
        ax.spines[s].set_visible(False)
    ax.tick_params(axis='y', labelsize=FS_TICK)

def draw_point_with_cis(ax, i):
    # 90% CI (solid)
    # ax.vlines(i, coef[i]-ci90[i], coef[i]+ci90[i], color=(0, 0, 0, ALPHA), lw=LW_90, zorder=2)
    # ax.hlines([coef[i]-ci90[i], coef[i]+ci90[i]], i-cap_hw, i+cap_hw, color=(0, 0, 0, ALPHA), lw=LW_90)
    # 95% CI (dashed)
    # ax.vlines(i, coef[i]-ci95[i], coef[i]-ci90[i], color=(0, 0, 0, ALPHA), ls='--', lw=LW_95)
    # ax.vlines(i, coef[i]+ci90[i], coef[i]+ci95[i], color=(0, 0, 0, ALPHA), ls='--', lw=LW_95)
    ax.vlines(i, coef[i]-ci95[i], coef[i]+ci95[i], color=(0, 0, 0, ALPHA), ls='--', lw=LW_95)

    ax.hlines([coef[i]-ci95[i], coef[i]+ci95[i]], i-cap_hw, i+cap_hw, color=(0, 0, 0, ALPHA), ls='--', lw=LW_95)
    # point & label
    ax.plot(i, coef[i], 'o', color='black', markersize=MS)
    ax.text(i+0.075, coef[i], f"{coef[i]:.2f}", ha='left', va='center', fontsize=FS_LAB)

def annotate_delta_line(ax, left_idx, right_idx, y, pair_slot):
    """Horizontal Δ line with caps and annotation."""
    d, se, st = pair_delta[pair_slot], pair_se[pair_slot], pair_star[pair_slot]
    # main line
    ax.plot([left_idx, right_idx], [y, y], 'k-', lw=1.8)
    # caps
    ax.vlines([left_idx, right_idx], y-0.03, y+0.03, color='k', lw=1.8)
    # annotation text (stars right after Δ value)
    xm = 0.5 * (left_idx + right_idx)
    ax.text(xm, y + 0.06, f"Δ = {d:+.2f}{st}  ({se:.2f})",
            ha='center', va='bottom', fontsize=FS_ANN)

def add_errorbar_legend(ax):
    """Add legend explaining 90% and 95% error bars — boxed, bottom-left."""
    from matplotlib.lines import Line2D
    legend_elems = [
        # Line2D([0], [0], color=(0, 0, 0, ALPHA), lw=LW_90, label='90% CI'),
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




def box_footer(ax, y0=-0.18, height=0.10):
    """
    Draw a boxed band under the x-axis in axes coords.
    y0, height are in axes coordinates (so y=0 is the x-axis).
    """
    rect = FancyBboxPatch(
        (0.0, y0),              # x,y (axes coords)
        1.0, height,            # width=full axis, height (axes coords)
        transform=ax.transAxes,
        boxstyle="round,pad=0.01",
        facecolor="0.97",       # light gray fill; use 'white' if you prefer
        edgecolor="0.6",
        linewidth=1.2,
        clip_on=False,          # <- important so it renders outside the axes
        zorder=0.3
    )
    ax.add_patch(rect)


# def add_bin_footer(ax, heading_text, values):
#     """
#     One-row footer with a light bounding box:
#       - Heading at axes origin (x=0 in axes coords)
#       - Numbers under tick marks
#       - Box encloses the whole row
#     """
#     y_row = -0.12           # vertical position in axes coordinates
#     box_height = 0.09       # height of the box (adjust as needed)
#     box_y = y_row - 0.02    # bottom edge of the box

#     # --- Background box ---
#     ax.add_patch(Rectangle(
#         (0, box_y),               # lower-left corner in axes coords
#         1, box_height,            # width=full axes width, height of box
#         transform=ax.transAxes,
#         facecolor='white', edgecolor='0.6',
#         lw=1.2, zorder=0, alpha=1.0
#     ))

#     # --- Heading anchored to left spine ---
#     ax.text(0.0, y_row, heading_text,
#             transform=ax.transAxes, ha='left', va='top',
#             fontsize=FS_FOOT, color='0.1', clip_on=False)

#     # --- Numbers aligned under bins ---
#     trans_blend = transforms.blended_transform_factory(ax.transData, ax.transAxes)
#     for xi, val in enumerate(values):
#         ax.text(xi, y_row, f"{val}",
#                 transform=trans_blend, ha='center', va='top',
#                 fontsize=FS_FOOT, color='0.1', clip_on=False)


def add_bin_footer(ax, heading_text, values, heading_x=-0.06, y_row=-0.12, visible_idx=None):
    """
    Footer with selective number visibility:
      - heading at axes x = heading_x (axes coords)
      - numbers under ticks (x=0,1,2 in data coords)
      - only show those indices in visible_idx
    """
    if visible_idx is None:
        visible_idx = range(len(values))

    # Heading anchored in AXES coords (so we can push it left of the spine)
    # ax.text(heading_x, y_row+.01, heading_text,
    #         transform=ax.transAxes, ha='left', va='top',
    #         fontsize=FS_FOOT, color='0.1', clip_on=False)
    ax.text(heading_x, y_row, heading_text,
            transform=ax.transAxes, ha='left', va='top',
            fontsize=FS_FOOT, color='0.1', clip_on=False)

    # Numbers: x in DATA coords (0,1,2), y in AXES coords (same y_row)
    trans_blend = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for xi, val in enumerate(values):
        if xi in visible_idx:
            ax.text(xi, y_row, f"{val}\%",
                    transform=trans_blend, ha='center', va='top',
                    fontsize=FS_FOOT, color='0.1', clip_on=False)



def set_visible_xlabels(ax, bins, visible_idx, fs):
    """Show labels only at indices in visible_idx; keep tick positions unchanged."""
    labels = [bins[i] if i in visible_idx else '' for i in range(len(bins))]
    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels(labels, fontsize=fs)


# ========== Data & display settings ==========

full_data_dict = {
    'BMI': {'coef': np.array([1.77, 1.05, 0.69]), #np.array([1.74, 1.05, 0.69]),
            'se': np.array([0.50, 0.58, 0.29]), #np.array([0.490, 0.58, 0.29]),
            'pair_delta': np.array([-0.72, -1.08, -0.36]), #np.array([-0.69, -1.05, -0.36]),
            'pair_se': np.array([0.59, 0.46, 0.54]), #np.array([0.59, 0.46, 0.54]),
            'pair_star': ['', '**', '']},
    'FIT': {'coef': np.array([3.92, 3.33, 2.69]),
            'se': np.array([0.52, 0.30, 0.26]),
            'pair_delta': np.array([-0.59, -1.23, -0.64]),
            'pair_se': np.array([0.35, 0.44, 0.24]),
            'pair_star': ['*', '***', '***']},
    'OFI': {'coef': np.array([3.13, 2.03, 0.91]), #np.array([4.67, 4.00, 2.78]),
            'se': np.array([0.10, 0.07, 0.06]), #np.array([0.17, 0.16, 0.16]),
            'pair_delta': np.array([-1.10, -2.22, -1.12]), #np.array([-0.66, -1.89, -1.22]),
            'pair_se': np.array([0.07, 0.08, 0.05]), #np.array([0.10, 0.13, 0.12]),
            'pair_star': ['***', '***', '***']},
}

pair_idx = [(0, 1), (0, 2), (1, 2)]
# bins = ['$[0,\sigma_t)$', '$[\sigma_t, 2 \sigma_t)$', '$[2\sigma_t, \infty)$']
bins = ['Small', 'Medium', 'Large']

# Figure aesthetics
FS_YLAB = 34
FS_TITLE = 40
FS_TICK = 30
FS_LAB  = 30
FS_ANN  = 30
FS_LEG  = 24
FS_FOOT = 28   # <- footer font size
MS = 16
LW_90 = 3.2
LW_95 = 1.8
ALPHA= 0.35
cap_hw = 0.03
x = np.arange(len(bins))


footer_values = {
    'BMI': [0.20, 1.08, 2.05],
    'FIT': [0.14, 0.63, 1.36],
    'OFI': [1.33, 5.67,13.65] #[0.6 , 2.17, 4.74],
}


# footer_heading = r"Avg. $|d_{n,t}|$:"
footer_heading = r"Avg. Size:"

names_dict = {
    'BMI': 'Benchmarking Intensity (BMI)',
    'FIT': 'Flow Induced Trading (FIT)',
    'OFI': 'Order Flow Imbalance (OFI)',
}

# ========== Loop through datasets ==========
for key, value in full_data_dict.items():
    this_dict = value
    coef = this_dict['coef']
    se = this_dict['se']
    pair_delta = this_dict['pair_delta']
    pair_se = this_dict['pair_se']
    pair_star = this_dict['pair_star']
    ci90 = 1.645 * se
    ci95 = 1.96  * se

    # --- Step 1 ---
    fig, ax = plt.subplots(figsize=(14, 10))
    base_plot(ax)
    set_visible_xlabels(ax, bins, visible_idx=[0], fs=FS_TICK)   # show only 'Small'
    draw_point_with_cis(ax, 0)
    add_errorbar_legend(ax)
    # box_footer(ax)
    add_bin_footer(ax, footer_heading, footer_values[key], visible_idx=[0])#, heading_x=-.2)
    # leave extra bottom space for footer
    plt.subplots_adjust(bottom=0.22)
    plt.title(names_dict[key], fontsize=FS_TITLE, pad=50)
    plt.savefig(f"step1_{key}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- Step 2 ---
    fig, ax = plt.subplots(figsize=(14, 10))
    base_plot(ax)
    set_visible_xlabels(ax, bins, visible_idx=[0, 1], fs=FS_TICK)   # show only 'Small' and 'Medium'
    for i in [0, 1]:
        draw_point_with_cis(ax, i)
    yline = max(coef[0]+ci95[0], coef[1]+ci95[1]) + 0.20
    annotate_delta_line(ax, *pair_idx[0], yline, pair_slot=0)
    add_errorbar_legend(ax)
    # box_footer(ax)
    add_bin_footer(ax, footer_heading, footer_values[key], visible_idx=[0, 1])#, heading_x=-.2)
    plt.subplots_adjust(bottom=0.22)
    plt.title(names_dict[key], fontsize=FS_TITLE, pad=50)
    plt.savefig(f"step2_{key}.png", dpi=300, bbox_inches='tight')
    plt.close()

    # --- Step 3 ---
    fig, ax = plt.subplots(figsize=(14, 10))
    base_plot(ax)
    set_visible_xlabels(ax, bins, visible_idx=[0,1,2], fs=FS_TICK)   # show only 'Small', 'Medium', and 'Large'

    for i in [0, 1, 2]:
        draw_point_with_cis(ax, i)
    y_big_small = max(coef[0]+ci95[0], coef[2]+ci95[2]) + 0.25
    y_big_medium = y_big_small - 0.65
    annotate_delta_line(ax, *pair_idx[1], y_big_small, pair_slot=1)
    annotate_delta_line(ax, *pair_idx[2], y_big_medium, pair_slot=2)
    add_errorbar_legend(ax)
    # box_footer(ax)
    add_bin_footer(ax, footer_heading, footer_values[key], visible_idx=[0, 1, 2])#, heading_x=-.2)
    plt.subplots_adjust(bottom=0.22)
    plt.title(names_dict[key], fontsize=FS_TITLE, pad=50)
    plt.savefig(f"step3_{key}.png", dpi=300, bbox_inches='tight')
    plt.close()
