import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ["PATH"] = "/usr/local/texlive/2020/bin/x86_64-darwin:" + os.environ["PATH"]

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pyreadr
from matplotlib.lines import Line2D


mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "text.latex.preamble": r"""
\usepackage[utf8]{inputenc}
\usepackage{sfmath}
\usepackage{newunicodechar}
\newunicodechar{Δ}{$\Delta$}
\newunicodechar{−}{-}
""",
})


SPECS = {
    "v3": {
        "path": "fm_stdev_posneg_v3.RDS",
        "title_suffix": "v3 pooled bins",
        "out_prefix": "asym_v3",
    },
    "v3_balanced": {
        "path": "fm_stdev_posneg_v3_balanced.RDS",
        "title_suffix": "v3 balanced months",
        "out_prefix": "asym_v3_balanced",
    },
    "v4": {
        "path": "fm_stdev_posneg_rebuilt_v4.RDS",
        "title_suffix": "v4 rebuilt bins",
        "out_prefix": "asym_v4",
    },
}

NAMES = {
    "BMI": "Benchmarking Intensity (BMI)",
    "FIT": "Flow Induced Trading (FIT)",
    "OFI": "Order Flow Imbalance (OFI)",
}

BINS = ["Medium - Small", "Big - Medium", "Big - Small"]
DIFF_VARS = {
    "pos": [
        "ofi_bin2_pos - ofi_bin1_pos",
        "ofi_bin3_pos - ofi_bin2_pos",
        "ofi_bin3_pos - ofi_bin1_pos",
    ],
    "neg": [
        "ofi_bin2_neg - ofi_bin1_neg",
        "ofi_bin3_neg - ofi_bin2_neg",
        "ofi_bin3_neg - ofi_bin1_neg",
    ],
}
BIN_VARS = {
    "pos": ["ofi_bin1_pos", "ofi_bin2_pos", "ofi_bin3_pos"],
    "neg": ["ofi_bin1_neg", "ofi_bin2_neg", "ofi_bin3_neg"],
}


FS_YLAB = 34
FS_TITLE = 38
FS_TICK = 30
FS_LAB = 26
FS_LEG = 24
FS_BOTTOM = 20
MS = 16
LW_95 = 1.8
ALPHA = 0.35
CAP_HW = 0.03
OFFSET = 0.20


def read_rds(path):
    result = pyreadr.read_r(path)
    return result[None]


def pick(df, type_name, var_name):
    rows = df[(df["type"] == type_name) & (df["var"] == var_name)]
    if len(rows) != 1:
        raise ValueError(f"Expected one row for {type_name=} {var_name=}; got {len(rows)}")
    return rows.iloc[0]


def base_plot(ax, y_max):
    ax.set_xlim(-0.5, 2.5)
    ax.set_ylim(-y_max * 1.6, y_max * 0.3)
    ax.set_xticks(range(3))
    ax.set_xticklabels(BINS, fontsize=FS_TICK)
    ax.set_ylabel(r"Multiplier Difference", fontsize=FS_YLAB)
    ax.axhline(0, color="0.7", linestyle="-", linewidth=1.5, zorder=1)
    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)
    ax.tick_params(axis="y", labelsize=FS_TICK)


def draw_point_with_cis(ax, i, coef_val, se_val, color, offset=0, bin_labels=None):
    ci95 = 1.96 * se_val
    x = i + offset
    ax.vlines(x, coef_val - ci95, coef_val + ci95,
              color=color, ls="--", lw=LW_95, alpha=ALPHA)
    ax.hlines([coef_val - ci95, coef_val + ci95],
              x - CAP_HW, x + CAP_HW,
              color=color, ls="--", lw=LW_95, alpha=ALPHA)
    ax.plot(x, coef_val, "o", color=color, markersize=MS)
    ax.text(x + 0.08, coef_val, f"{coef_val:.1f}",
            ha="left", va="center", fontsize=FS_LAB, color=color)
    if bin_labels is not None:
        ax.text(x, coef_val - ci95 - 0.2, bin_labels[i],
                ha="center", va="top", fontsize=FS_BOTTOM, color=color)


def add_errorbar_legend(ax):
    legend_elems = [
        Line2D([0], [0], color="black", marker="o", linestyle="None",
               markersize=MS, label="Positive"),
        Line2D([0], [0], color="red", marker="o", linestyle="None",
               markersize=MS, label="Negative"),
        Line2D([0], [0], color=(0, 0, 0, ALPHA), lw=LW_95, ls="--", label=r"95\% CI"),
    ]
    leg = ax.legend(handles=legend_elems,
                    fontsize=FS_LEG,
                    loc="lower left",
                    frameon=True,
                    facecolor="white",
                    edgecolor="0.6",
                    framealpha=1.0,
                    handlelength=2.8,
                    bbox_to_anchor=(0.02, 0.02))
    for text in leg.get_texts():
        text.set_color("0.1")


def labels_from_bins(bin_coefs):
    return [
        f"{bin_coefs[1]:.1f} - {bin_coefs[0]:.1f}",
        f"{bin_coefs[2]:.1f} - {bin_coefs[1]:.1f}",
        f"{bin_coefs[2]:.1f} - {bin_coefs[0]:.1f}",
    ]


def plot_one_spec(spec_name, spec):
    df = read_rds(spec["path"])
    df = df[df["var_added"] == "controls_char+controls_liq"].copy()
    out_dir = os.path.join("plots", spec_name)
    os.makedirs(out_dir, exist_ok=True)

    for type_name in ["OFI", "FIT", "BMI"]:
        pos_bins = [pick(df, type_name, var_name)["coef"] for var_name in BIN_VARS["pos"]]
        neg_bins = [pick(df, type_name, var_name)["coef"] for var_name in BIN_VARS["neg"]]
        pos_labels = labels_from_bins(pos_bins)
        neg_labels = labels_from_bins(neg_bins)

        pos_rows = [pick(df, type_name, var_name) for var_name in DIFF_VARS["pos"]]
        neg_rows = [pick(df, type_name, var_name) for var_name in DIFF_VARS["neg"]]
        pos_coefs = np.array([row["coef"] for row in pos_rows])
        pos_ses = np.array([row["se"] for row in pos_rows])
        neg_coefs = np.array([row["coef"] for row in neg_rows])
        neg_ses = np.array([row["se"] for row in neg_rows])

        all_vals = np.concatenate([
            pos_coefs - 1.96 * pos_ses,
            pos_coefs + 1.96 * pos_ses,
            neg_coefs - 1.96 * neg_ses,
            neg_coefs + 1.96 * neg_ses,
        ])
        y_max = max(np.max(np.abs(all_vals)), 0.1)

        fig, ax = plt.subplots(figsize=(14, 10))
        base_plot(ax, y_max)
        for i in range(3):
            draw_point_with_cis(ax, i, pos_coefs[i], pos_ses[i], "black", -OFFSET, pos_labels)
            draw_point_with_cis(ax, i, neg_coefs[i], neg_ses[i], "red", OFFSET, neg_labels)

        add_errorbar_legend(ax)
        plt.subplots_adjust(bottom=0.15)
        title = f"{NAMES[type_name]}: {spec['title_suffix']}"
        plt.title(title, fontsize=FS_TITLE, pad=50)
        out_path = os.path.join(out_dir, f"{spec['out_prefix']}_{type_name}.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(out_path)


if __name__ == "__main__":
    for spec_name, spec in SPECS.items():
        plot_one_spec(spec_name, spec)
