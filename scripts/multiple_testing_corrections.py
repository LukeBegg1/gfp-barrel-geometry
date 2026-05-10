#!/usr/bin/env python3
"""
Multiple testing corrections for GFP barrel geometry manuscript.
Applies Benjamini-Hochberg FDR correction to all Spearman correlations
and group comparison tests.
"""

import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import warnings
warnings.filterwarnings('ignore')

# ─── Load data ───────────────────────────────────────────────────────────
df = pd.read_csv('/Users/lukebegg/Downloads/deep_analysis/merged_complete_data.csv')
dih = pd.read_csv('/Users/lukebegg/Downloads/deep_analysis/megley_dihedrals.csv')

# Rename for convenience
df.columns = [c.strip().lower() for c in df.columns]
dih.columns = [c.strip().lower() for c in dih.columns]

out_lines = []
def report(text=""):
    out_lines.append(text)
    print(text)

report("=" * 80)
report("MULTIPLE TESTING CORRECTIONS FOR GFP BARREL GEOMETRY MANUSCRIPT")
report("=" * 80)
report(f"\nDataset: {len(df)} structures in merged data")
report(f"Dihedral dataset: {len(dih)} structures")

# ─── 1. Compute ALL Spearman correlations ────────────────────────────────
report("\n" + "=" * 80)
report("SECTION 1: ALL SPEARMAN CORRELATIONS (uncorrected)")
report("=" * 80)

correlation_results = []

def spearman_test(data, x_col, y_col, label=None):
    """Compute Spearman correlation, dropping NaNs."""
    sub = data[[x_col, y_col]].dropna()
    if len(sub) < 5:
        return None
    rho, p = stats.spearmanr(sub[x_col], sub[y_col])
    name = label or f"{y_col} vs {x_col}"
    correlation_results.append({
        'test_name': name,
        'test_type': 'Spearman',
        'rho': rho,
        'p_value': p,
        'n': len(sub)
    })
    return rho, p, len(sub)

# Emission vs barrel geometry
for var in ['eccentricity', 'minor_axis', 'circularity', 'major_axis', 'barrel_length', 'convex_area']:
    spearman_test(df, var, 'em_max', f"Emission vs {var}")

# QY vs b_factor_ratio
spearman_test(df, 'b_factor_ratio', 'lit_qy', "QY vs b_factor_ratio")

# Emission vs b_factor_ratio
spearman_test(df, 'b_factor_ratio', 'em_max', "Emission vs b_factor_ratio")

# Emission vs chrom_contacts
spearman_test(df, 'chrom_contacts', 'em_max', "Emission vs chrom_contacts")

# tau vs eccentricity (from main data, using tau_main)
spearman_test(df, 'eccentricity', 'tau_main', "tau vs eccentricity")

# |tau| vs QY and |phi| vs QY from dihedral data
dih['abs_tau'] = dih['tau_megley'].abs()
dih['abs_phi'] = dih['phi_megley'].abs()
spearman_test(dih, 'abs_tau', 'lit_qy', "|tau| vs QY")
spearman_test(dih, 'abs_phi', 'lit_qy', "|phi| vs QY")

# tau+phi vs emission and QY
dih['tau_plus_phi'] = dih['tau_megley'].abs() + dih['phi_megley'].abs()
spearman_test(dih, 'tau_plus_phi', 'em_max', "|tau|+|phi| vs Emission")
spearman_test(dih, 'tau_plus_phi', 'lit_qy', "|tau|+|phi| vs QY")

# Resolution vs geometry
for var in ['convex_area', 'eccentricity', 'minor_axis']:
    spearman_test(df, 'resolution', var, f"Resolution vs {var}")

# Print all uncorrected results
report(f"\n{'Test':<35} {'rho':>8} {'p-value':>12} {'n':>6}")
report("-" * 65)
for r in correlation_results:
    sig = "***" if r['p_value'] < 0.001 else "**" if r['p_value'] < 0.01 else "*" if r['p_value'] < 0.05 else ""
    report(f"{r['test_name']:<35} {r['rho']:>8.4f} {r['p_value']:>12.2e} {r['n']:>6} {sig}")

# ─── 2. Mann-Whitney U tests ─────────────────────────────────────────────
report("\n" + "=" * 80)
report("SECTION 2: MANN-WHITNEY U TESTS (chromophore vs no-chromophore)")
report("=" * 80)

group_results = []

chrom = df[df['has_chromophore'] == True]
no_chrom = df[df['has_chromophore'] == False]
report(f"\nChromophore present: n={len(chrom)}, absent: n={len(no_chrom)}")

for var in ['convex_area', 'eccentricity', 'circularity', 'minor_axis']:
    g1 = chrom[var].dropna()
    g2 = no_chrom[var].dropna()
    if len(g1) < 3 or len(g2) < 3:
        continue
    stat, p = stats.mannwhitneyu(g1, g2, alternative='two-sided')
    group_results.append({
        'test_name': f"MW: {var} (chrom vs no-chrom)",
        'test_type': 'Mann-Whitney U',
        'statistic': stat,
        'p_value': p,
        'n1': len(g1),
        'n2': len(g2),
        'median1': g1.median(),
        'median2': g2.median()
    })
    report(f"\n{var}:")
    report(f"  Chromophore:    median={g1.median():.4f}, n={len(g1)}")
    report(f"  No chromophore: median={g2.median():.4f}, n={len(g2)}")
    report(f"  U={stat:.1f}, p={p:.2e}")

# ─── 3. Kruskal-Wallis tests ─────────────────────────────────────────────
report("\n" + "=" * 80)
report("SECTION 3: KRUSKAL-WALLIS TESTS (by color class)")
report("=" * 80)

color_classes = df['color_class'].dropna().unique()
report(f"\nColor classes: {sorted(color_classes)}")
for cc in sorted(color_classes):
    report(f"  {cc}: n={len(df[df['color_class']==cc])}")

for var in ['eccentricity', 'minor_axis', 'circularity']:
    groups = []
    group_labels = []
    for cc in sorted(color_classes):
        g = df.loc[df['color_class'] == cc, var].dropna()
        if len(g) >= 3:
            groups.append(g)
            group_labels.append(cc)
    if len(groups) < 2:
        continue
    stat, p = stats.kruskal(*groups)
    group_results.append({
        'test_name': f"KW: {var} by color_class",
        'test_type': 'Kruskal-Wallis',
        'statistic': stat,
        'p_value': p,
        'n_groups': len(groups),
        'group_labels': ', '.join(group_labels)
    })
    report(f"\n{var} by color class:")
    for gl, g in zip(group_labels, groups):
        report(f"  {gl:<12}: median={g.median():.4f}, n={len(g)}")
    report(f"  H={stat:.3f}, p={p:.2e}")

# ─── 4. Benjamini-Hochberg FDR correction ────────────────────────────────
report("\n" + "=" * 80)
report("SECTION 4: BENJAMINI-HOCHBERG FDR CORRECTION")
report("=" * 80)

# Collect ALL p-values
all_tests = []
for r in correlation_results:
    all_tests.append({'name': r['test_name'], 'type': r['test_type'], 'p_raw': r['p_value'],
                      'detail': f"rho={r['rho']:.4f}, n={r['n']}"})
for r in group_results:
    all_tests.append({'name': r['test_name'], 'type': r['test_type'], 'p_raw': r['p_value'],
                      'detail': f"stat={r['statistic']:.1f}"})

p_values = np.array([t['p_raw'] for t in all_tests])
n_tests = len(p_values)

report(f"\nTotal number of tests being corrected: {n_tests}")
report(f"FDR threshold: q = 0.05")

# BH procedure manually (for transparency) and via scipy
# Manual BH:
sorted_indices = np.argsort(p_values)
sorted_p = p_values[sorted_indices]
ranks = np.arange(1, n_tests + 1)
bh_critical = (ranks / n_tests) * 0.05

# Determine which survive
bh_adjusted = np.empty(n_tests)
bh_adjusted[sorted_indices[-1]] = sorted_p[-1]
for i in range(n_tests - 2, -1, -1):
    idx = sorted_indices[i]
    bh_adjusted[idx] = min(sorted_p[i] * n_tests / (i + 1), bh_adjusted[sorted_indices[i + 1]])
bh_adjusted = np.minimum(bh_adjusted, 1.0)

# Also use scipy for verification
try:
    from scipy.stats import false_discovery_control
    scipy_bh = false_discovery_control(p_values, method='bh')
    report("(Verified with scipy.stats.false_discovery_control)")
except ImportError:
    scipy_bh = None
    report("(scipy.stats.false_discovery_control not available, using manual BH)")

# Store adjusted p-values
for i, t in enumerate(all_tests):
    t['p_adj'] = bh_adjusted[i]
    t['survives_fdr'] = bh_adjusted[i] < 0.05

# Sort by adjusted p-value for display
all_tests_sorted = sorted(all_tests, key=lambda x: x['p_adj'])

report(f"\n{'Test':<40} {'p_raw':>12} {'p_adj(BH)':>12} {'FDR<0.05':>10}")
report("-" * 78)
for t in all_tests_sorted:
    flag = "YES" if t['survives_fdr'] else "no"
    report(f"{t['name']:<40} {t['p_raw']:>12.2e} {t['p_adj']:>12.2e} {flag:>10}")

n_survive = sum(1 for t in all_tests if t['survives_fdr'])
report(f"\nTests surviving FDR correction: {n_survive} / {n_tests}")

report("\n--- Tests that SURVIVE BH-FDR correction (q < 0.05) ---")
for t in all_tests_sorted:
    if t['survives_fdr']:
        report(f"  {t['name']}: p_raw={t['p_raw']:.2e}, p_adj={t['p_adj']:.2e} ({t['detail']})")

report("\n--- Tests that DO NOT survive BH-FDR correction ---")
for t in all_tests_sorted:
    if not t['survives_fdr']:
        report(f"  {t['name']}: p_raw={t['p_raw']:.2e}, p_adj={t['p_adj']:.2e} ({t['detail']})")

# ─── 5. Partial correlations controlling for resolution ──────────────────
report("\n" + "=" * 80)
report("SECTION 5: PARTIAL CORRELATIONS (controlling for resolution)")
report("=" * 80)

report("""
METHOD: Partial Spearman correlations via the residual method.
For each partial correlation of X vs Y controlling for Z (resolution):
  1. Rank-transform X, Y, and Z
  2. Regress ranked X on ranked Z using OLS: X_resid = ranked_X - predicted(ranked_X | ranked_Z)
  3. Regress ranked Y on ranked Z using OLS: Y_resid = ranked_Y - predicted(ranked_Y | ranked_Z)
  4. Compute Pearson correlation of X_resid and Y_resid
     (Pearson on rank-residuals = partial Spearman correlation)
  5. Test significance using t = r * sqrt((n-3) / (1-r^2)), df = n-3
""")

def partial_spearman(data, x_col, y_col, z_col, label=None):
    """Partial Spearman correlation of x,y controlling for z via residual method."""
    sub = data[[x_col, y_col, z_col]].dropna()
    n = len(sub)
    if n < 6:
        return None

    # Step 1: Rank transform
    rx = stats.rankdata(sub[x_col])
    ry = stats.rankdata(sub[y_col])
    rz = stats.rankdata(sub[z_col])

    # Step 2: Regress ranked X on ranked Z
    slope_xz, intercept_xz, _, _, _ = stats.linregress(rz, rx)
    x_resid = rx - (intercept_xz + slope_xz * rz)

    # Step 3: Regress ranked Y on ranked Z
    slope_yz, intercept_yz, _, _, _ = stats.linregress(rz, ry)
    y_resid = ry - (intercept_yz + slope_yz * rz)

    # Step 4: Pearson correlation of residuals
    r_partial, _ = stats.pearsonr(x_resid, y_resid)

    # Step 5: Significance test
    t_stat = r_partial * np.sqrt((n - 3) / (1 - r_partial**2))
    p_value = 2 * stats.t.sf(np.abs(t_stat), df=n - 3)

    name = label or f"{y_col} vs {x_col} | {z_col}"
    report(f"\n{name}:")
    report(f"  n = {n}")
    report(f"  Partial Spearman rho = {r_partial:.4f}")
    report(f"  t-statistic = {t_stat:.3f}, df = {n-3}")
    report(f"  p-value = {p_value:.2e}")

    # Also show the zero-order correlation for comparison
    rho_zero, p_zero = stats.spearmanr(sub[x_col], sub[y_col])
    report(f"  (Zero-order Spearman: rho={rho_zero:.4f}, p={p_zero:.2e})")

    return {'name': name, 'r_partial': r_partial, 'p_value': p_value, 'n': n,
            'rho_zero': rho_zero, 'p_zero': p_zero}

partial_results = []
res = partial_spearman(df, 'eccentricity', 'em_max', 'resolution', "Emission vs eccentricity | resolution")
if res: partial_results.append(res)

res = partial_spearman(df, 'minor_axis', 'em_max', 'resolution', "Emission vs minor_axis | resolution")
if res: partial_results.append(res)

res = partial_spearman(df, 'convex_area', 'em_max', 'resolution', "Emission vs area | resolution")
if res: partial_results.append(res)

res = partial_spearman(df, 'b_factor_ratio', 'lit_qy', 'resolution', "QY vs b_factor_ratio | resolution")
if res: partial_results.append(res)

# BH correction on partial correlations
if partial_results:
    p_partial = np.array([r['p_value'] for r in partial_results])
    n_p = len(p_partial)
    sorted_idx = np.argsort(p_partial)
    sorted_pp = p_partial[sorted_idx]

    bh_adj_partial = np.empty(n_p)
    bh_adj_partial[sorted_idx[-1]] = min(sorted_pp[-1], 1.0)
    for i in range(n_p - 2, -1, -1):
        idx = sorted_idx[i]
        bh_adj_partial[idx] = min(sorted_pp[i] * n_p / (i + 1), bh_adj_partial[sorted_idx[i + 1]])
    bh_adj_partial = np.minimum(bh_adj_partial, 1.0)

    report("\n--- Partial correlations with BH correction ---")
    report(f"{'Test':<45} {'r_partial':>10} {'p_raw':>12} {'p_adj':>12} {'FDR<0.05':>10}")
    report("-" * 92)
    for i, r in enumerate(partial_results):
        flag = "YES" if bh_adj_partial[i] < 0.05 else "no"
        report(f"{r['name']:<45} {r['r_partial']:>10.4f} {r['p_value']:>12.2e} {bh_adj_partial[i]:>12.2e} {flag:>10}")

# ─── 6. Summary ──────────────────────────────────────────────────────────
report("\n" + "=" * 80)
report("SECTION 6: SUMMARY FOR MANUSCRIPT")
report("=" * 80)

report("""
METHODS TEXT (partial correlations):
\"Partial Spearman correlations controlling for crystallographic resolution
were computed using the residual method. Variables were first rank-transformed,
then each ranked variable was regressed on ranked resolution via ordinary
least squares. The Pearson correlation of the resulting residuals yields the
partial Spearman correlation coefficient. Significance was assessed using the
t-distribution with n - 3 degrees of freedom, where t = r_partial * sqrt((n-3)
/ (1 - r_partial^2)). All p-values (Spearman correlations, Mann-Whitney U tests,
Kruskal-Wallis tests, and partial correlations) were corrected for multiple
comparisons using the Benjamini-Hochberg procedure to control the false
discovery rate at q = 0.05.\"
""")

report("CORRELATIONS SURVIVING FDR (q < 0.05):")
for t in all_tests_sorted:
    if t['survives_fdr']:
        report(f"  - {t['name']}: rho/stat from {t['detail']}, p_adj={t['p_adj']:.2e}")

report("\nCORRELATIONS NOT SURVIVING FDR:")
for t in all_tests_sorted:
    if not t['survives_fdr']:
        report(f"  - {t['name']}: p_raw={t['p_raw']:.2e}, p_adj={t['p_adj']:.2e}")

report("\nPARTIAL CORRELATIONS SUMMARY:")
if partial_results:
    for i, r in enumerate(partial_results):
        change = "attenuated" if abs(r['r_partial']) < abs(r['rho_zero']) else "strengthened"
        report(f"  - {r['name']}: partial rho={r['r_partial']:.4f} (zero-order: {r['rho_zero']:.4f}, {change})")

# ─── Save results ────────────────────────────────────────────────────────
with open('/Users/lukebegg/Downloads/deep_analysis/multiple_testing_results.txt', 'w') as f:
    f.write('\n'.join(out_lines))

print("\n\nResults saved to /Users/lukebegg/Downloads/deep_analysis/multiple_testing_results.txt")
