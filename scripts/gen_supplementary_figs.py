#!/usr/bin/env python3
"""Generate supplementary Figures S1-S5 for JCIM GFP barrel geometry manuscript."""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# ── Global style ──────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'axes.linewidth': 1.0,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'savefig.facecolor': 'white',
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# ── Color scheme ──────────────────────────────────────────────────────────────
CLASS_COLORS = {
    'cyan':   '#00ACC1',
    'green':  '#43A047',
    'yellow': '#FDD835',
    'orange': '#FB8C00',
    'red':    '#E53935',
    'blue':   '#1E88E5',
}

OUTDIR = '/Users/lukebegg/Downloads/deep_analysis/pub_figures'

# ── Load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv('/Users/lukebegg/Downloads/deep_analysis/merged_complete_data.csv')
print("Columns:", list(df.columns))
print(f"Loaded {len(df)} rows")


# ── Helper functions ──────────────────────────────────────────────────────────
def clean_spines(ax):
    """Remove top and right spines."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def panel_label(ax, label, x=-0.12, y=1.08):
    """Add bold panel label."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=14, fontweight='bold', va='top', ha='left')

def safe_ylim(ax, data, pad_frac=0.05):
    """Set ylim robustly, handling NaN/Inf."""
    clean = data.replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) == 0:
        return
    lo, hi = clean.min(), clean.max()
    pad = (hi - lo) * pad_frac if hi != lo else 1.0
    ax.set_ylim(lo - pad, hi + pad)

def safe_xlim(ax, data, pad_frac=0.05):
    clean = data.replace([np.inf, -np.inf], np.nan).dropna()
    if len(clean) == 0:
        return
    lo, hi = clean.min(), clean.max()
    pad = (hi - lo) * pad_frac if hi != lo else 1.0
    ax.set_xlim(lo - pad, hi + pad)

def spearman_text(x, y):
    """Return formatted Spearman rho and p."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5:
        return 'n < 5'
    rho, p = stats.spearmanr(x[mask], y[mask])
    if p < 0.001:
        return f'ρ = {rho:.2f}, p < 0.001'
    else:
        return f'ρ = {rho:.2f}, p = {p:.3f}'


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S1: Emission vs Geometry  (1×3)
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(7.5, 6))

pairs = [
    ('minor_axis',   'Minor Axis (Å)'),
    ('eccentricity', 'Eccentricity'),
    ('circularity',  'Circularity'),
]
labels = ['(A)', '(B)', '(C)']

sub = df.dropna(subset=['em_max'])

for ax, (col, xlabel), lab in zip(axes, pairs, labels):
    for cls in ['cyan', 'green', 'yellow', 'orange', 'red', 'blue']:
        sel = sub[sub['color_class'] == cls]
        sel_valid = sel.dropna(subset=[col])
        if len(sel_valid) == 0:
            continue
        ax.scatter(sel_valid[col], sel_valid['em_max'],
                   c=CLASS_COLORS[cls], s=20, alpha=0.6,
                   edgecolors='none', label=cls.capitalize(), zorder=2)
    clean_spines(ax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Emission Max (nm)')
    txt = spearman_text(sub[col].values, sub['em_max'].values)
    ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=9,
            va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='grey', alpha=0.8, pad=3))
    panel_label(ax, lab)
    safe_ylim(ax, sub['em_max'])

# Single legend below
handles, lbls = axes[0].get_legend_handles_labels()
fig.legend(handles, lbls, loc='lower center', ncol=6, frameon=False,
           fontsize=9, bbox_to_anchor=(0.5, -0.02))

fig.tight_layout(rect=[0, 0.05, 1, 1])
path = f'{OUTDIR}/fig3_emission_vs_geometry.png'
fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {path}')


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S2: Cis/Trans Configuration  (2×2)
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(7.5, 6))
axes = axes.flatten()

cfg_sub = df[df['config'].isin(['cis', 'trans', 'twisted'])].copy()
cfg_order = ['cis', 'trans', 'twisted']
cfg_colors = {'cis': '#43A047', 'trans': '#FB8C00', 'twisted': '#E53935'}

metrics = [
    ('eccentricity', 'Eccentricity'),
    ('circularity',  'Circularity'),
    ('minor_axis',   'Minor Axis (Å)'),
    ('convex_area',  'Area (Å²)'),
]
labels_s2 = ['(A)', '(B)', '(C)', '(D)']

for ax, (col, ylabel), lab in zip(axes, metrics, labels_s2):
    data_groups = [cfg_sub.loc[cfg_sub['config'] == c, col].dropna() for c in cfg_order]
    bp = ax.boxplot(data_groups, positions=[1, 2, 3], widths=0.6,
                    patch_artist=True, showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.4))
    for patch, c in zip(bp['boxes'], cfg_order):
        patch.set_facecolor(cfg_colors[c])
        patch.set_alpha(0.7)
    for element in ['whiskers', 'caps', 'medians']:
        for item in bp[element]:
            item.set_color('black')
            item.set_linewidth(0.8)
    for item in bp['medians']:
        item.set_linewidth(1.5)

    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['Cis', 'Trans', 'Twisted'])
    ax.set_ylabel(ylabel)
    clean_spines(ax)
    panel_label(ax, lab)
    safe_ylim(ax, cfg_sub[col])

    # Mann-Whitney cis vs trans
    cis_d = data_groups[0]
    tra_d = data_groups[1]
    if len(cis_d) >= 3 and len(tra_d) >= 3:
        stat, p = stats.mannwhitneyu(cis_d, tra_d, alternative='two-sided')
        if p < 0.001:
            ptxt = 'p < 0.001'
        else:
            ptxt = f'p = {p:.3f}'
        ymax = max(cis_d.max(), tra_d.max())
        yrange = ax.get_ylim()[1] - ax.get_ylim()[0]
        bar_y = ymax + 0.05 * yrange
        ax.plot([1, 1, 2, 2], [bar_y, bar_y + 0.02*yrange,
                                bar_y + 0.02*yrange, bar_y],
                color='black', linewidth=0.8)
        ax.text(1.5, bar_y + 0.03*yrange, ptxt, ha='center', va='bottom', fontsize=8)

fig.tight_layout()
path = f'{OUTDIR}/fig9_cis_trans.png'
fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {path}')


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S3: Resolution Control  (1×3)
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(7.5, 6))

res_sub = df.dropna(subset=['resolution', 'minor_axis']).copy()
res_sub = res_sub[np.isfinite(res_sub['resolution']) & np.isfinite(res_sub['minor_axis'])]

# (A) Resolution vs Minor Axis
ax = axes[0]
ax.scatter(res_sub['resolution'], res_sub['minor_axis'],
           c='#1E88E5', s=20, alpha=0.5, edgecolors='none')
clean_spines(ax)
ax.set_xlabel('Resolution (Å)')
ax.set_ylabel('Minor Axis (Å)')
txt = spearman_text(res_sub['resolution'].values, res_sub['minor_axis'].values)
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(facecolor='white', edgecolor='grey', alpha=0.8, pad=3))
panel_label(ax, '(A)')

# (B) Emission vs Minor Axis colored by resolution
ax = axes[1]
em_sub = res_sub.dropna(subset=['em_max'])
em_sub = em_sub[np.isfinite(em_sub['em_max'])]
sc = ax.scatter(em_sub['minor_axis'], em_sub['em_max'],
                c=em_sub['resolution'], cmap='viridis', s=20, alpha=0.6,
                edgecolors='none')
cbar = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label('Resolution (Å)', fontsize=9)
cbar.ax.tick_params(labelsize=8)
clean_spines(ax)
ax.set_xlabel('Minor Axis (Å)')
ax.set_ylabel('Emission Max (nm)')
txt = spearman_text(em_sub['minor_axis'].values, em_sub['em_max'].values)
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(facecolor='white', edgecolor='grey', alpha=0.8, pad=3))
panel_label(ax, '(B)')

# (C) High-resolution subset (<2.0 Å)
ax = axes[2]
hires = em_sub[em_sub['resolution'] < 2.0]
ax.scatter(hires['minor_axis'], hires['em_max'],
           c='#43A047', s=20, alpha=0.6, edgecolors='none')
clean_spines(ax)
ax.set_xlabel('Minor Axis (Å)')
ax.set_ylabel('Emission Max (nm)')
txt = spearman_text(hires['minor_axis'].values, hires['em_max'].values)
ax.text(0.05, 0.95, txt, transform=ax.transAxes, fontsize=9, va='top',
        bbox=dict(facecolor='white', edgecolor='grey', alpha=0.8, pad=3))
ax.text(0.05, 0.85, f'n = {len(hires)} (< 2.0 Å)', transform=ax.transAxes,
        fontsize=9, va='top',
        bbox=dict(facecolor='white', edgecolor='grey', alpha=0.8, pad=3))
panel_label(ax, '(C)')

fig.tight_layout()
path = f'{OUTDIR}/fig6_resolution.png'
fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {path}')


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S4: Correlation Heatmap  (square)
# ══════════════════════════════════════════════════════════════════════════════
hm_cols = ['convex_area', 'minor_axis', 'major_axis', 'eccentricity',
           'circularity', 'em_max', 'lit_qy', 'b_factor_ratio',
           'chrom_contacts', 'resolution']
hm_labels = ['Area', 'Minor', 'Major', 'Ecc.', 'Circ.',
             'Em.', 'QY', 'B-ratio', 'Contacts', 'Resol.']

hm_data = df[hm_cols].copy()
n = len(hm_cols)
rho_mat = np.full((n, n), np.nan)
p_mat = np.full((n, n), np.nan)

for i in range(n):
    for j in range(n):
        mask = np.isfinite(hm_data.iloc[:, i]) & np.isfinite(hm_data.iloc[:, j])
        if mask.sum() >= 5:
            r, p = stats.spearmanr(hm_data.iloc[:, i][mask], hm_data.iloc[:, j][mask])
            rho_mat[i, j] = r
            p_mat[i, j] = p

fig, ax = plt.subplots(figsize=(8, 8))
im = ax.imshow(rho_mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='equal')
cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label('Spearman ρ', fontsize=11)

for i in range(n):
    for j in range(n):
        val = rho_mat[i, j]
        pv = p_mat[i, j]
        if np.isnan(val):
            continue
        stars = ''
        if pv < 0.001:
            stars = '***'
        elif pv < 0.01:
            stars = '**'
        elif pv < 0.05:
            stars = '*'
        color = 'white' if abs(val) > 0.5 else 'black'
        ax.text(j, i, f'{val:.2f}{stars}', ha='center', va='center',
                fontsize=9, color=color)

ax.set_xticks(range(n))
ax.set_yticks(range(n))
ax.set_xticklabels(hm_labels, rotation=45, ha='right', fontsize=10)
ax.set_yticklabels(hm_labels, fontsize=10)
ax.tick_params(top=False, bottom=True, left=True, right=False)

fig.tight_layout()
path = f'{OUTDIR}/fig7_heatmap.png'
fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {path}')


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE S5: Chromophore Types  (3×1 stacked)
# ══════════════════════════════════════════════════════════════════════════════
ct_counts = df['chromophore_type'].value_counts()
valid_types = ct_counts[ct_counts >= 5].index.tolist()
ct_sub = df[df['chromophore_type'].isin(valid_types)].copy()

# Sort by median minor axis
medians = ct_sub.groupby('chromophore_type')['minor_axis'].median().sort_values()
sorted_types = medians.index.tolist()

# Labels with n
type_labels = []
for t in sorted_types:
    nn = ct_counts[t]
    type_labels.append(f'{t} (n={nn})')

fig, axes = plt.subplots(3, 1, figsize=(7.5, 8))

metrics_s5 = [
    ('minor_axis',   'Minor Axis (Å)'),
    ('eccentricity', 'Eccentricity'),
    ('em_max',       'Emission Max (nm)'),
]
labels_s5 = ['(A)', '(B)', '(C)']

box_color = '#1E88E5'

for ax, (col, ylabel), lab in zip(axes, metrics_s5, labels_s5):
    data_groups = [ct_sub.loc[ct_sub['chromophore_type'] == t, col].dropna()
                   for t in sorted_types]
    # Filter out empty groups
    non_empty = [(g, l) for g, l in zip(data_groups, type_labels) if len(g) > 0]
    if len(non_empty) == 0:
        continue
    groups_clean = [g for g, l in non_empty]
    labels_clean = [l for g, l in non_empty]

    bp = ax.boxplot(groups_clean, positions=range(1, len(groups_clean)+1),
                    widths=0.6, patch_artist=True, showfliers=True,
                    flierprops=dict(marker='o', markersize=3, alpha=0.4))
    for patch in bp['boxes']:
        patch.set_facecolor(box_color)
        patch.set_alpha(0.6)
    for element in ['whiskers', 'caps', 'medians']:
        for item in bp[element]:
            item.set_color('black')
            item.set_linewidth(0.8)
    for item in bp['medians']:
        item.set_linewidth(1.5)

    ax.set_xticks(range(1, len(labels_clean)+1))
    ax.set_xticklabels(labels_clean, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel(ylabel)
    clean_spines(ax)
    panel_label(ax, lab)
    all_vals = pd.concat(groups_clean)
    safe_ylim(ax, all_vals)

fig.tight_layout()
path = f'{OUTDIR}/fig8_chromophore_types.png'
fig.savefig(path, dpi=300, bbox_inches='tight', facecolor='white')
plt.close(fig)
print(f'Saved: {path}')

print('\nAll 5 supplementary figures generated successfully.')
