#!/usr/bin/env python3
"""
SENSITIVITY ANALYSIS for GFP barrel cross-section geometry.
Tests robustness to: slice thickness (1.5, 2.5 vs 2.0 A) and backbone-only atoms.
"""

import numpy as np
import csv
from pathlib import Path
from scipy.spatial import ConvexHull
from scipy.linalg import eigh
from scipy.spatial.transform import Rotation
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

import gemmi

STRUCT_DIR = Path("/Users/lukebegg/Downloads/scop_gfp_structures")
MERGED_CSV = Path("/Users/lukebegg/Downloads/deep_analysis/merged_complete_data.csv")
ORIGINAL_CSV = Path("/Users/lukebegg/Downloads/cross_section_correct.csv")
OUTPUT_FILE = Path("/Users/lukebegg/Downloads/deep_analysis/sensitivity_results.txt")

CHROMOPHORE_RESIDUES = {'CRO', 'CR2', 'GYS', 'SYG', 'CRQ', 'CRF', 'CRW', 'CRY',
                        'NRQ', 'NYG', 'CH6', 'CH7', 'CRG', 'CRU', 'CRV', 'CRS',
                        '66A', 'CR0', 'GYC', 'LYG', 'TYG', 'OHD', 'SWG', 'QYG',
                        'CR7', 'CR8', 'CR9', 'CRK', 'RC7', 'CFY', 'PIA', 'B2H'}

STANDARD_AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
               'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'MSE', 'PHE', 'PRO',
               'SER', 'THR', 'TRP', 'TYR', 'VAL'}

BACKBONE_ATOMS = {'N', 'CA', 'C', 'O'}


def rotation_matrix_from_vectors(vec1, vec2):
    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)
    if np.allclose(a, b):
        return np.eye(3)
    if np.allclose(a, -b):
        ortho = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        ortho = ortho - np.dot(ortho, a) * a
        ortho = ortho / np.linalg.norm(ortho)
        return Rotation.from_rotvec(np.pi * ortho).as_matrix()
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))
    return R


def analyze_structure_variant(cif_file, pdb_id, thickness=2.0, backbone_only=False):
    """Analyze structure with given slice thickness and atom selection."""
    try:
        structure = gemmi.read_structure(str(cif_file))
        model = structure[0]

        best_chain = None
        best_count = 0
        for chain in model:
            count = sum(1 for res in chain if res.name in STANDARD_AA)
            if count > best_count:
                best_count = count
                best_chain = chain

        if best_chain is None or best_count < 100:
            return None

        ca_atoms = []
        all_atoms = []
        chromophore_atoms = []
        chromophore_type = None

        for residue in best_chain:
            for atom in residue:
                if atom.element.name == 'H':
                    continue
                # For backbone_only mode, only keep N, CA, C, O from standard residues
                if backbone_only and residue.name in STANDARD_AA:
                    if atom.name not in BACKBONE_ATOMS:
                        continue

                pos = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                all_atoms.append(pos)

                if residue.name in STANDARD_AA and atom.name == 'CA':
                    ca_atoms.append(pos)

            if residue.name in CHROMOPHORE_RESIDUES:
                for atom in residue:
                    if atom.element.name != 'H':
                        chromophore_atoms.append(np.array([atom.pos.x, atom.pos.y, atom.pos.z]))
                chromophore_type = residue.name

        ca_atoms = np.array(ca_atoms)
        all_atoms = np.array(all_atoms)

        if len(ca_atoms) < 100:
            return None

        # PCA on CA atoms
        ca_centroid = np.mean(ca_atoms, axis=0)
        ca_centered = ca_atoms - ca_centroid
        cov = np.cov(ca_centered.T)
        eigenvalues, eigenvectors = eigh(cov)
        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        barrel_axis = eigenvectors[:, 0]
        if barrel_axis[2] < 0:
            barrel_axis = -barrel_axis

        # Rotate to z-axis
        z_axis = np.array([0, 0, 1])
        R = rotation_matrix_from_vectors(barrel_axis, z_axis)
        all_atoms_rotated = np.dot(all_atoms - ca_centroid, R.T)

        if chromophore_atoms:
            chromophore_atoms = np.array(chromophore_atoms)
            chromophore_rotated = np.dot(chromophore_atoms - ca_centroid, R.T)
            center_z = np.mean(chromophore_rotated, axis=0)[2]
        else:
            center_z = np.mean(all_atoms_rotated[:, 2])

        all_atoms_final = all_atoms_rotated.copy()
        all_atoms_final[:, 2] -= center_z

        # Slice
        mask = np.abs(all_atoms_final[:, 2]) <= thickness
        slice_atoms = all_atoms_final[mask]

        if len(slice_atoms) < 10:
            return None

        points_2d = slice_atoms[:, :2]

        try:
            hull = ConvexHull(points_2d)
            area = hull.volume
            perimeter = hull.area
        except Exception:
            return None

        cov_2d = np.cov(points_2d.T)
        eigvals_2d, _ = eigh(cov_2d)
        eigvals_2d = np.sort(eigvals_2d)[::-1]

        major_axis = 4 * np.sqrt(eigvals_2d[0])
        minor_axis = 4 * np.sqrt(eigvals_2d[1])
        eccentricity = np.sqrt(1 - (minor_axis / major_axis) ** 2) if major_axis > 0 else 0
        circularity = 4 * np.pi * area / (perimeter ** 2) if perimeter > 0 else 0

        return {
            'pdb_id': pdb_id,
            'convex_area': area,
            'major_axis': major_axis,
            'minor_axis': minor_axis,
            'eccentricity': eccentricity,
            'circularity': circularity,
            'n_atoms': len(slice_atoms),
        }

    except Exception as e:
        return None


def select_representative_structures():
    """Select ~50 representative structures spanning all color classes."""
    with open(MERGED_CSV) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    # Group by color class
    by_color = {}
    for r in rows:
        cc = r['color_class']
        if cc not in by_color:
            by_color[cc] = []
        by_color[cc].append(r)

    print("Color class distribution:")
    for cc, items in sorted(by_color.items(), key=lambda x: -len(x[1])):
        print(f"  {cc}: {len(items)}")

    # Target counts
    targets = {'green': 12, 'red': 10, 'yellow': 10, 'cyan': 10, 'orange': 5, '': 5}
    selected = []

    for cc, target in targets.items():
        pool = by_color.get(cc, [])
        # Prefer structures with spectral data
        with_spec = [r for r in pool if r['em_max'] and r['em_max'] != '']
        without_spec = [r for r in pool if not r['em_max'] or r['em_max'] == '']

        chosen = []
        # Take as many with spectral data as possible
        if len(with_spec) >= target:
            # Spread across the range by sorting by em_max and taking evenly spaced
            with_spec.sort(key=lambda r: float(r['em_max']))
            step = max(1, len(with_spec) // target)
            chosen = [with_spec[i] for i in range(0, len(with_spec), step)][:target]
        else:
            chosen = with_spec[:]
            remaining = target - len(chosen)
            chosen += without_spec[:remaining]

        selected.extend(chosen)

    pdb_ids = [r['pdb_id'] for r in selected]
    print(f"\nSelected {len(pdb_ids)} representative structures")
    return pdb_ids, selected


def main():
    print("=" * 70)
    print("GFP BARREL CROSS-SECTION SENSITIVITY ANALYSIS")
    print("=" * 70)

    # 1. Select representative structures
    pdb_ids, selected_rows = select_representative_structures()

    # Build lookup for original values and spectral data
    orig_lookup = {}
    with open(ORIGINAL_CSV) as f:
        reader = csv.DictReader(f)
        for r in reader:
            orig_lookup[r['pdb_id']] = r

    spectral_lookup = {}
    for r in selected_rows:
        spectral_lookup[r['pdb_id']] = r

    # 2. Run variants
    variants = {
        'thin_slice_1.5A': {'thickness': 1.5, 'backbone_only': False},
        'thick_slice_2.5A': {'thickness': 2.5, 'backbone_only': False},
        'backbone_only_2.0A': {'thickness': 2.0, 'backbone_only': True},
    }

    variant_results = {name: {} for name in variants}

    for i, pdb_id in enumerate(pdb_ids):
        cif_file = STRUCT_DIR / f"{pdb_id.lower()}.cif"
        if not cif_file.exists():
            print(f"  Missing CIF: {pdb_id}")
            continue

        for vname, vparams in variants.items():
            result = analyze_structure_variant(cif_file, pdb_id, **vparams)
            if result:
                variant_results[vname][pdb_id] = result

        if (i + 1) % 10 == 0:
            print(f"  Processed {i+1}/{len(pdb_ids)}...")

    print(f"\nProcessing complete.")

    # 3. Compare each variant to original
    metrics = ['convex_area', 'eccentricity', 'circularity', 'minor_axis', 'major_axis']
    output_lines = []

    def log(msg=""):
        print(msg)
        output_lines.append(msg)

    log("=" * 70)
    log("SENSITIVITY ANALYSIS RESULTS")
    log("=" * 70)
    log(f"Representative sample: {len(pdb_ids)} structures")
    log("")

    for vname, vdata in variant_results.items():
        log(f"\n{'=' * 70}")
        log(f"VARIANT: {vname}")
        log(f"  Structures successfully computed: {len(vdata)}")
        log(f"{'=' * 70}")

        # Match with original data
        common_ids = [pid for pid in vdata if pid in orig_lookup]
        log(f"  Structures with original comparison: {len(common_ids)}")

        if len(common_ids) < 5:
            log("  Too few common structures for comparison.")
            continue

        log(f"\n  {'Metric':<15} {'Pearson r':>10} {'Spearman rho':>13} {'MAD':>10} {'Mean Orig':>10} {'Mean Var':>10} {'% Change':>9}")
        log(f"  {'-'*78}")

        for metric in metrics:
            orig_vals = []
            var_vals = []
            for pid in common_ids:
                ov = orig_lookup[pid].get(metric, '')
                vv = vdata[pid].get(metric, '')
                if ov and vv and ov != '' and vv != '':
                    try:
                        orig_vals.append(float(ov))
                        var_vals.append(float(vv))
                    except (ValueError, TypeError):
                        pass

            if len(orig_vals) < 5:
                log(f"  {metric:<15} insufficient data")
                continue

            orig_arr = np.array(orig_vals)
            var_arr = np.array(var_vals)

            pearson_r, pearson_p = stats.pearsonr(orig_arr, var_arr)
            spearman_r, spearman_p = stats.spearmanr(orig_arr, var_arr)
            mad = np.mean(np.abs(orig_arr - var_arr))
            mean_orig = np.mean(orig_arr)
            mean_var = np.mean(var_arr)
            pct_change = (mean_var - mean_orig) / mean_orig * 100 if mean_orig != 0 else 0

            log(f"  {metric:<15} {pearson_r:>10.4f} {spearman_r:>13.4f} {mad:>10.3f} {mean_orig:>10.3f} {mean_var:>10.3f} {pct_change:>+8.2f}%")

    # 4. Emission correlations for each variant
    log(f"\n\n{'=' * 70}")
    log("EMISSION WAVELENGTH CORRELATIONS BY VARIANT")
    log("=" * 70)

    # Original correlations (from merged data)
    log("\n  ORIGINAL (2.0A all-atom):")
    em_ecc_orig = []
    em_minor_orig = []
    for r in selected_rows:
        em = r.get('em_max', '')
        ecc = r.get('eccentricity', '')
        mi = r.get('minor_axis', '')
        if em and em != '' and ecc and ecc != '':
            try:
                em_ecc_orig.append((float(em), float(ecc)))
            except (ValueError, TypeError):
                pass
        if em and em != '' and mi and mi != '':
            try:
                em_minor_orig.append((float(em), float(mi)))
            except (ValueError, TypeError):
                pass

    if len(em_ecc_orig) >= 5:
        em_vals, ecc_vals = zip(*em_ecc_orig)
        r, p = stats.pearsonr(em_vals, ecc_vals)
        rho, rho_p = stats.spearmanr(em_vals, ecc_vals)
        log(f"    em_max vs eccentricity:  Pearson r={r:.4f} (p={p:.4f}), Spearman rho={rho:.4f} (p={rho_p:.4f})  [n={len(em_ecc_orig)}]")

    if len(em_minor_orig) >= 5:
        em_vals, mi_vals = zip(*em_minor_orig)
        r, p = stats.pearsonr(em_vals, mi_vals)
        rho, rho_p = stats.spearmanr(em_vals, mi_vals)
        log(f"    em_max vs minor_axis:    Pearson r={r:.4f} (p={p:.4f}), Spearman rho={rho:.4f} (p={rho_p:.4f})  [n={len(em_minor_orig)}]")

    for vname, vdata in variant_results.items():
        log(f"\n  VARIANT: {vname}")
        em_ecc = []
        em_minor = []
        for r in selected_rows:
            pid = r['pdb_id']
            em = r.get('em_max', '')
            if not em or em == '' or pid not in vdata:
                continue
            try:
                em_f = float(em)
            except (ValueError, TypeError):
                continue
            vr = vdata[pid]
            em_ecc.append((em_f, vr['eccentricity']))
            em_minor.append((em_f, vr['minor_axis']))

        if len(em_ecc) >= 5:
            em_vals, ecc_vals = zip(*em_ecc)
            r, p = stats.pearsonr(em_vals, ecc_vals)
            rho, rho_p = stats.spearmanr(em_vals, ecc_vals)
            log(f"    em_max vs eccentricity:  Pearson r={r:.4f} (p={p:.4f}), Spearman rho={rho:.4f} (p={rho_p:.4f})  [n={len(em_ecc)}]")

        if len(em_minor) >= 5:
            em_vals, mi_vals = zip(*em_minor)
            r, p = stats.pearsonr(em_vals, mi_vals)
            rho, rho_p = stats.spearmanr(em_vals, mi_vals)
            log(f"    em_max vs minor_axis:    Pearson r={r:.4f} (p={p:.4f}), Spearman rho={rho:.4f} (p={rho_p:.4f})  [n={len(em_minor)}]")

    # 5. Summary assessment
    log(f"\n\n{'=' * 70}")
    log("ROBUSTNESS SUMMARY")
    log("=" * 70)
    log("")
    log("Interpretation guide:")
    log("  Pearson r > 0.95:  Excellent agreement (values are stable)")
    log("  Pearson r > 0.90:  Good agreement (minor sensitivity)")
    log("  Pearson r > 0.80:  Moderate agreement (some sensitivity)")
    log("  Pearson r < 0.80:  Poor agreement (parameter-sensitive)")
    log("")
    log("  Spearman rho > 0.95: Rank order strongly preserved")
    log("  Spearman rho > 0.90: Rank order well preserved")
    log("")

    # Save
    with open(OUTPUT_FILE, 'w') as f:
        f.write("\n".join(output_lines))
    print(f"\nResults saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
