#!/usr/bin/env python3
"""
VALIDATION SCRIPT — "second Claude code to check-up on itself"

Independently re-calculates cross-section geometry for the 17 wild-type FP
crystal structures from raw CIF files and compares to the stored values in
merged_complete_data.csv. Flags any discrepancy > tolerance.

This is the quality-control counterpart to the main analysis pipeline:
it uses a slightly different implementation (no shared code) to verify
the stored results are internally consistent.

Usage: python3 validate_analysis.py
"""

import numpy as np
import gemmi
import pandas as pd
from pathlib import Path
from scipy.spatial import ConvexHull

CIF_DIR = Path.home() / 'Downloads' / 'scop_gfp_structures'
MERGED  = Path.home() / 'Downloads' / 'deep_analysis' / 'merged_complete_data.csv'
OUT     = Path.home() / 'Downloads' / 'deep_analysis' / 'validation_report.txt'

# 17 wild-type structures to validate
WT_17 = ['1W7S','1GGX','2A46','1XAE','1XMZ',
         '1UIS','3PIB','3PJB','4EDO',
         '2GW3','1ZUX','2DDC',
         '4OHS','5LTQ','4JF9','5EXB','8I4J']

STANDARD_AA = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL'
}
CHROMOPHORE_RESIDUES = {
    'CRO','CRQ','CR7','CRK','CH7','CH6','CR2','CR8','QYG','GYS',
    'NRQ','5SQ','IIC','HTB','PYG','PYK','KYC','MYK','YCM','ASJ',
    'QMM','SEP','FGL','4IN','4TE','SYG','TYG','SRG'
}

SLICE_HALF  = 2.0
TOLERANCE   = 5.0    # Å or area units — flag if difference > this

lines = []
def log(s=''):
    print(s)
    lines.append(s)

log('=' * 70)
log('VALIDATION REPORT — FP barrel cross-section geometry')
log('Independent re-calculation vs stored values in merged_complete_data.csv')
log('=' * 70)
log()

# Load stored values
df = pd.read_csv(MERGED)
df['pdb_id'] = df['pdb_id'].str.upper()

def recompute(cif_path, pdb_id):
    """
    Independent cross-section calculation.
    Uses the SAME physical approach (PCA + convex hull + 4*sqrt covariance)
    but implemented from scratch to avoid shared-code bias.
    """
    try:
        st = gemmi.read_structure(str(cif_path))
        model = st[0]

        # Find chromophore chain/residue
        chrom_pos = None
        ca_all, heavy_all = [], []
        chain_residue_counts = {}

        for chain in model:
            n_aa = sum(1 for r in chain if r.name in STANDARD_AA)
            chain_residue_counts[chain.name] = n_aa

        best_chain_name = max(chain_residue_counts, key=chain_residue_counts.get)
        best_chain = model[best_chain_name]

        for res in best_chain:
            is_std = res.name in STANDARD_AA
            is_chrom = res.name in CHROMOPHORE_RESIDUES

            if is_std or is_chrom:
                for atom in res:
                    if not atom.is_hydrogen():
                        heavy_all.append([atom.pos.x, atom.pos.y, atom.pos.z])
            if is_std:
                ca = res.find_atom('CA', '\0')
                if ca:
                    ca_all.append([ca.pos.x, ca.pos.y, ca.pos.z])
            if is_chrom and chrom_pos is None:
                # Use centroid of chromophore heavy atoms
                chrom_atoms = []
                for atom in res:
                    if not atom.is_hydrogen():
                        chrom_atoms.append([atom.pos.x, atom.pos.y, atom.pos.z])
                if chrom_atoms:
                    chrom_pos = np.mean(chrom_atoms, axis=0)

        if len(ca_all) < 80:
            return None, "too few CA atoms"

        ca_arr    = np.array(ca_all)
        heavy_arr = np.array(heavy_all)

        # PCA on CA to get barrel axis
        centroid = ca_arr.mean(axis=0)
        cov3d    = np.cov((ca_arr - centroid).T)
        eigvals, eigvecs = np.linalg.eigh(cov3d)
        barrel_axis = eigvecs[:, np.argmax(eigvals)]

        # Barrel length
        ca_proj = (ca_arr - centroid) @ barrel_axis
        barrel_length = float(ca_proj.max() - ca_proj.min())

        # Rotation matrix: barrel_axis → Z
        zhat = np.array([0.0, 0.0, 1.0])
        v    = np.cross(barrel_axis, zhat)
        s    = np.linalg.norm(v)
        c    = np.dot(barrel_axis, zhat)
        if s < 1e-10:
            R = np.eye(3)
        else:
            vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
            R  = np.eye(3) + vx + vx @ vx * (1 - c) / s**2

        heavy_rot = (heavy_arr - centroid) @ R.T

        # Chromophore z-level
        if chrom_pos is not None:
            z_chrom = float(((chrom_pos - centroid) @ R.T)[2])
        else:
            z_chrom = 0.0   # fallback: barrel centre

        # Slice
        mask  = np.abs(heavy_rot[:, 2] - z_chrom) <= SLICE_HALF
        pts2d = heavy_rot[mask, :2]
        if len(pts2d) < 20:
            return None, f"only {len(pts2d)} atoms in slice"

        hull  = ConvexHull(pts2d)
        area  = float(hull.volume)
        perim = float(hull.area)
        circ  = 4 * np.pi * area / perim**2

        # Axes: 4*sqrt of covariance eigenvalues of ALL slice atoms
        cov2d  = np.cov(pts2d.T)
        ev2d,_ = np.linalg.eigh(cov2d)
        ev2d   = np.sort(np.abs(ev2d))[::-1]
        major  = 4 * np.sqrt(ev2d[0])
        minor  = 4 * np.sqrt(ev2d[1])
        ecc    = float(np.sqrt(1 - ev2d[1]/ev2d[0])) if ev2d[0] > 0 else 0.0

        return {
            'area':          round(area,   1),
            'major_axis':    round(major,  2),
            'minor_axis':    round(minor,  2),
            'eccentricity':  round(ecc,    3),
            'circularity':   round(circ,   3),
            'barrel_length': round(barrel_length, 1),
            'z_chrom':       round(z_chrom, 2),
            'n_slice':       int(mask.sum()),
        }, None

    except Exception as e:
        return None, str(e)


# ── Run validation ────────────────────────────────────────────────────────────
pass_count = fail_count = warn_count = 0

log(f"{'PDB':6} {'Metric':14} {'Stored':>10} {'Recomputed':>10} {'Diff':>8} {'Status':>8}")
log('-' * 62)

for pdb_id in WT_17:
    cif_path = CIF_DIR / f'{pdb_id.lower()}.cif'
    if not cif_path.exists():
        log(f'{pdb_id:6}  CIF not found at {cif_path}')
        fail_count += 1
        continue

    recomp, err = recompute(cif_path, pdb_id)
    if recomp is None:
        log(f'{pdb_id:6}  RECOMPUTE FAILED: {err}')
        fail_count += 1
        continue

    stored_row = df[df['pdb_id'] == pdb_id]
    if stored_row.empty:
        log(f'{pdb_id:6}  NOT IN merged_complete_data.csv')
        fail_count += 1
        continue

    stored = stored_row.iloc[0]
    checks = [
        ('area',         'convex_area',  TOLERANCE * 10),   # area in Å²
        ('major_axis',   'major_axis',   TOLERANCE),
        ('minor_axis',   'minor_axis',   TOLERANCE),
        ('eccentricity', 'eccentricity', 0.05),
        ('circularity',  'circularity',  0.05),
        ('barrel_length','barrel_length',TOLERANCE),
    ]

    pdb_passed = True
    for re_col, st_col, tol in checks:
        rv = recomp[re_col]
        sv = stored.get(st_col, np.nan)
        if pd.isna(sv):
            status = 'STORED=NaN'
        else:
            diff = abs(rv - sv)
            if diff > tol:
                status = '⚠ MISMATCH'
                warn_count += 1
                pdb_passed = False
            else:
                status = '✓ OK'
                pass_count += 1
        log(f'{pdb_id:6}  {re_col:14} {sv:>10.3f}  {rv:>10.3f}  {abs(rv-sv) if not pd.isna(sv) else 0:>8.3f}  {status}')

    log()

log('=' * 62)
log(f'Checks passed: {pass_count}')
log(f'Warnings:      {warn_count}')
log(f'Failures:      {fail_count}')
log()

if warn_count == 0 and fail_count == 0:
    log('RESULT: ALL CHECKS PASSED — stored values are consistent with')
    log('        independent re-calculation from raw CIF files.')
elif warn_count > 0:
    log('RESULT: SOME MISMATCHES DETECTED — review flagged values above.')
    log('        Differences may reflect rounding, alt conformers, or')
    log('        implementation differences in the ellipse-fitting step.')
else:
    log('RESULT: STRUCTURAL FAILURES — check CIF file paths.')

OUT.write_text('\n'.join(lines))
log(f'\nReport saved to: {OUT}')
