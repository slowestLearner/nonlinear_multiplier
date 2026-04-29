import pyreadr
import pandas as pd
import numpy as np
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))
INPUT_DIR = os.path.join('..', '..', '20250117_quarterly')

# ========== Step 1: Load and describe ==========

pos_df = pyreadr.read_r(os.path.join(INPUT_DIR, 'fm_stdev_pos.RDS'))[None]
neg_df = pyreadr.read_r(os.path.join(INPUT_DIR, 'fm_stdev_neg.RDS'))[None]

print("pos shape:", pos_df.shape)
print("neg shape:", neg_df.shape)
print("types:", sorted(pos_df['type'].unique()))
print("spec_idx:", sorted(pos_df['spec_idx'].unique()))
print("vars (sample):", list(pos_df['var'].unique())[:8])

# ========== Step 2: Filter ==========

DIFF_VARS = ['ofi_bin2 - ofi_bin1', 'ofi_bin3 - ofi_bin2', 'ofi_bin3 - ofi_bin1']
SPEC_IDX = [1, 2, 3]

def filter_data(df, label):
    out = df[(df['spec_idx'].isin(SPEC_IDX)) & (df['var'].isin(DIFF_VARS))].copy()
    print(f"{label}: {len(out)} rows after filter (expect 27 = 3 types × 3 specs × 3 diffs)")
    return out

pos_filt = filter_data(pos_df, 'pos')
neg_filt = filter_data(neg_df, 'neg')

# ========== Step 3: Validate ==========

for label, df in [('pos', pos_filt), ('neg', neg_filt)]:
    expected = {(t, s, v) for t in ['BMI', 'FIT', 'OFI']
                for s in SPEC_IDX for v in DIFF_VARS}
    actual = set(zip(df['type'], df['spec_idx'], df['var']))
    missing = expected - actual
    print(f"{label}: {len(missing)} missing combos: {missing}")
    print(f"  NaN coef: {df['coef'].isna().sum()}, NaN se: {df['se'].isna().sum()}")
    ofi3 = df[(df['type'] == 'OFI') & (df['spec_idx'] == 3)][['var', 'coef', 'se']]
    print(f"  OFI spec3:\n{ofi3.to_string(index=False)}")

# ========== Step 4: Format and build table ==========

THRESHOLDS = [2.576, 1.960, 1.645]

def format_coef(coef, se):
    abs_t = abs(coef / se)
    n_stars = sum(abs_t > t for t in THRESHOLDS)
    stars = '^{' + '*' * n_stars + '}$' if n_stars > 0 else '$'
    return f'${coef:.2f}{stars}'

def format_se(se):
    return f'$({se:.2f})$'

ROW_LABELS = {
    'ofi_bin2 - ofi_bin1': (r'$M_{\{|d_{n,t}| \in [\sigma, 2\sigma]\}}'
                             r' - M_{\{|d_{n,t}| < \sigma\}}$'),
    'ofi_bin3 - ofi_bin2': (r'$M_{\{|d_{n,t}| > 2 \sigma\}}'
                             r' - M_{\{|d_{n,t}| \in [\sigma, 2\sigma]\}}$'),
    'ofi_bin3 - ofi_bin1': (r'$M_{\{|d_{n,t}| > 2 \sigma\}}'
                             r' - M_{\{|d_{n,t}| < \sigma\}}$'),
}
TYPES = ['BMI', 'FIT', 'OFI']


def build_panel_rows(df):
    rows = []
    for var in DIFF_VARS:
        coef_cells, se_cells = [], []
        for type_name in TYPES:
            for spec in SPEC_IDX:
                r = df[(df['type'] == type_name) & (df['spec_idx'] == spec) & (df['var'] == var)]
                c, s = r['coef'].values[0], r['se'].values[0]
                coef_cells.append(format_coef(c, s))
                se_cells.append(format_se(s))
        rows.append((ROW_LABELS[var], coef_cells, se_cells))
    return rows


def fmt_row(label, cells, vspace=False):
    v = r' \vspace{5pt}' if vspace else ''
    return f'  {label} & {" & ".join(cells)} \\\\{v} \n'


def build_tabular(pos_rows, neg_rows):
    lines = [r'\begin{tabular}{lccccccccc}' + '\n']

    for panel_label, rows in [('Panel A: Positive shocks', pos_rows),
                               ('Panel B: Negative shocks', neg_rows)]:
        lines.append(f'  \\hline \\multicolumn{{10}}{{c}}{{{panel_label}}} \\\\\n')
        lines.append(' \\hline\n')
        lines.append(r'  & \multicolumn{3}{c}{BMI} & \multicolumn{3}{c}{FIT} & \multicolumn{3}{c}{OFI} \\' + ' \n')
        lines.append(r'  \cmidrule(l){2-4} \cmidrule(l){5-7} \cmidrule(l){8-10}' + ' \n')
        lines.append('  & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) \\\\ \n')
        for i, (label, coef_cells, se_cells) in enumerate(rows):
            is_last = (i == len(rows) - 1)
            lines.append(fmt_row(label, coef_cells, vspace=is_last))
            lines.append(fmt_row('', se_cells))

    lines.append('   \\hline \n')
    lines.append(r'\end{tabular}' + '\n')
    return ''.join(lines)


pos_rows = build_panel_rows(pos_filt)
neg_rows = build_panel_rows(neg_filt)
tabular = build_tabular(pos_rows, neg_rows)

# ========== Step 5: Write output ==========

out_path = 'reg_asym_stdev.tex'
with open(out_path, 'w') as f:
    f.write(tabular)

print(f"\nWritten to {out_path} ({len(tabular)} chars)")
print("\n--- Output ---")
print(tabular)
