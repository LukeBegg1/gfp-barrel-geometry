#!/usr/bin/env python3
"""Multivariate analysis for GFP barrel geometry manuscript."""

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
import statsmodels.api as sm
import sys, os

out_path = "/Users/lukebegg/Downloads/deep_analysis/multivariate_results.txt"
lines = []
def p(s=""):
    print(s)
    lines.append(s)

# ── Load data ──
main = pd.read_csv("/Users/lukebegg/Downloads/deep_analysis/merged_complete_data.csv")
dihed = pd.read_csv("/Users/lukebegg/Downloads/deep_analysis/megley_dihedrals.csv")

p("=== COLUMN NAMES ===")
p(f"Main data columns ({len(main.columns)}): {list(main.columns)}")
p(f"Dihedral columns ({len(dihed.columns)}): {list(dihed.columns)}")
p(f"Main data rows: {len(main)}")
p(f"Dihedral rows: {len(dihed)}")
p()

# ── Merge ──
# Dihedral file has many columns that overlap; only take pdb_id and tau_megley
dihed_sub = dihed[['pdb_id', 'tau_megley']].copy()
df = main.merge(dihed_sub, on='pdb_id', how='inner')
p(f"After merge: {len(df)} rows")

# Compute |tau|
df['abs_tau'] = df['tau_megley'].abs()

# ── Helper: run OLS ──
def run_ols(df, outcome_col, predictor_cols, label):
    p(f"\n{'='*60}")
    p(f"  MULTIPLE REGRESSION: {label}")
    p(f"  Outcome: {outcome_col}")
    p(f"  Predictors: {predictor_cols}")
    p(f"{'='*60}")

    cols_needed = predictor_cols + [outcome_col]
    sub = df[cols_needed].dropna()
    p(f"  N (complete cases): {len(sub)}")

    if len(sub) < 10:
        p("  *** Too few cases for regression ***")
        return

    y = sub[outcome_col].values
    X_raw = sub[predictor_cols].values

    # Standardize predictors
    scaler = StandardScaler()
    X_std = scaler.fit_transform(X_raw)

    # Add constant for intercept
    X_std_c = sm.add_constant(X_std)

    model = sm.OLS(y, X_std_c).fit()

    p(f"\n  R²          = {model.rsquared:.4f}")
    p(f"  Adjusted R² = {model.rsquared_adj:.4f}")
    p(f"  F-statistic = {model.fvalue:.2f}, p = {model.f_pvalue:.2e}")
    p(f"  N = {int(model.nobs)}, df_resid = {int(model.df_resid)}")
    p()

    p(f"  {'Predictor':<20s} {'Std Beta':>10s} {'SE':>10s} {'t':>10s} {'p-value':>12s} {'Sig':>5s}")
    p(f"  {'-'*67}")
    for i, name in enumerate(predictor_cols):
        idx = i + 1  # skip constant
        beta = model.params[idx]
        se = model.bse[idx]
        t = model.tvalues[idx]
        pv = model.pvalues[idx]
        sig = "***" if pv < 0.001 else "**" if pv < 0.01 else "*" if pv < 0.05 else ""
        p(f"  {name:<20s} {beta:>10.4f} {se:>10.4f} {t:>10.3f} {pv:>12.4e} {sig:>5s}")

    p(f"\n  Intercept: {model.params[0]:.4f}")

    # VIF check
    p(f"\n  Correlation matrix of predictors:")
    corr = pd.DataFrame(X_raw, columns=predictor_cols).corr()
    for c in predictor_cols:
        vals = " ".join([f"{corr.loc[c,c2]:>6.2f}" for c2 in predictor_cols])
        p(f"    {c:<20s} {vals}")

    return model


# ── Helper: Random Forest ──
def run_rf(df, outcome_col, predictor_cols, label):
    p(f"\n{'='*60}")
    p(f"  RANDOM FOREST: {label}")
    p(f"  Outcome: {outcome_col}")
    p(f"{'='*60}")

    cols_needed = predictor_cols + [outcome_col]
    sub = df[cols_needed].dropna()
    p(f"  N (complete cases): {len(sub)}")

    if len(sub) < 10:
        p("  *** Too few cases ***")
        return

    y = sub[outcome_col].values
    X = sub[predictor_cols].values

    rf = RandomForestRegressor(n_estimators=500, oob_score=True, random_state=42, n_jobs=-1)
    rf.fit(X, y)

    p(f"\n  OOB R² = {rf.oob_score_:.4f}")

    # MDI importances
    p(f"\n  Feature importances (Mean Decrease Impurity):")
    p(f"  {'Predictor':<20s} {'Importance':>12s}")
    p(f"  {'-'*34}")
    mdi = rf.feature_importances_
    order = np.argsort(mdi)[::-1]
    for i in order:
        p(f"  {predictor_cols[i]:<20s} {mdi[i]:>12.4f}")

    # Permutation importances
    perm = permutation_importance(rf, X, y, n_repeats=30, random_state=42, n_jobs=-1)
    p(f"\n  Permutation importances (30 repeats):")
    p(f"  {'Predictor':<20s} {'Mean':>10s} {'Std':>10s}")
    p(f"  {'-'*42}")
    order2 = np.argsort(perm.importances_mean)[::-1]
    for i in order2:
        p(f"  {predictor_cols[i]:<20s} {perm.importances_mean[i]:>10.4f} {perm.importances_std[i]:>10.4f}")


# ── Run analyses ──
predictors = ['b_factor_ratio', 'eccentricity', 'minor_axis', 'circularity', 'abs_tau', 'resolution']

# 1. QY regression
run_ols(df, 'lit_qy', predictors, "Predicting Quantum Yield (lit_qy)")

# 2. Emission regression
run_ols(df, 'em_max', predictors, "Predicting Emission Maximum (em_max)")

# 3. Random forest for QY
run_rf(df, 'lit_qy', predictors, "Predicting Quantum Yield (lit_qy)")

# ── Save ──
with open(out_path, 'w') as f:
    f.write("\n".join(lines))
p(f"\nResults saved to {out_path}")
