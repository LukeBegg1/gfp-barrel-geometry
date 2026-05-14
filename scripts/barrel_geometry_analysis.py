#!/usr/bin/env python3
"""
CORRECT CROSS-SECTION ANALYSIS
==============================
Properly aligned perpendicular to barrel axis.

Method:
1. Calculate barrel axis (PCA of Cα atoms)
2. Rotate structure so barrel axis = z-axis
3. Center on chromophore (or barrel center)
4. Slice at z=0 ± 2Å (now perpendicular to barrel)

This gives TRUE cross-sections, not oblique cuts.
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import ConvexHull
from scipy.linalg import eigh
from scipy.spatial.transform import Rotation
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

try:
    import gemmi
except ImportError:
    print("ERROR: pip3 install gemmi")
    exit(1)

REMOVE_IDS = {'9N6A', '7T9X', '6JC5', '9YZQ', '6B0B', '6BBO', '8XPP', '7BWN'}
NO_CHROMOPHORE = {
    # Manually verified structures lacking a mature chromophore
    # (reduced set after correcting residue detection; ~41 structures total)
    '1QXT', '1QY3', '1QYO', '1S6Z',
    '1YHG', '1YHH', '1YHI', '1Z1P', '1Z1Q', '2AWJ',
    '2G2S', '2G3D', '2G5Z',
    '2HFC', '4GES', '4GF6', '4JFG', '5B61', '6M9Y',
    '7SF9', '8AAB', '8C0N'
}

# Imidazolinone scaffold atom names (fallback detection)
SCAFFOLD_ATOMS = {'CA2', 'CB2', 'CG2', 'N2', 'C2', 'O2', 'N3', 'C3'}

# Known non-chromophore modified residues
NON_CHROMOPHORE_MODS = {'CSD', 'DHA', 'E1H', 'LEN', 'NFA', 'NLW', 'NME', 'PQ4',
                        'MSE', 'SEP', 'TPO', 'PTR', 'CME', 'OCS'}

CHROMOPHORE_RESIDUES = {
    # Validated: imidazolinone ring + aromatic ring at pos 66 (Tyr/Trp/His-derived)
    # Tyr66-based (GFP, YFP, RFP, and variants)
    'CRO', 'CR2', 'CRQ', 'CRF', 'CRW', 'CRY', 'CRG', 'CRU', 'CRV', 'CRS',
    '66A', 'CR0', 'GYC', 'LYG', 'TYG', 'OHD', 'SWG', 'QYG', 'NRQ', 'NYG',
    'CH6', 'CH7', 'RC7', 'B2H', 'GYG', 'GYS',
    'CR7', 'CR9', 'CRK', 'CRH', 'CRI', 'CRJ', 'CRL', 'CRM', 'CRN', 'CRP',
    'CRT', 'CRZ', 'CRB', 'CRC', 'CRD', 'CRE',
    'XYG', '5SQ', 'C12', 'IO8', 'KXV', 'KY7', 'KY4', 'KZ1', 'KZY',
    'BJO', 'BJF', 'QLG', 'KZ7', 'KZ4', 'KZG', 'KZV', 'OFM', '4M9',
    'IEY', 'DYG', '7R0', 'QIP', 'QYX', 'JBY', 'BF6', 'CIV',
    'CJO', 'CQ1', 'CQ2', 'CZO', 'EYG', 'FHE', 'GMO',
    'LKE', 'M3V', 'MCQ', 'MFC', 'OIM', 'QCA', 'TUK', 'VUB', 'X9Q',
    'XXY', '0YG', '7R6', 'HMF', 'MYG', 'RCQ', 'YGS', 'ZYG',
    # Trp66-based (CFP-type)
    'PIA', 'CFY', 'CCY',
    # His66-based (BFP-type)
    'IIC', 'CSH',
    # EXCLUDED (cyclized but no aromatic at pos 66, not fluorescent):
    # CR8 (Arg66 mutagenesis), CRX, MDO, CWR, C99, CLV, KWS, VYA, Q2K,
    # NRP, CR5 (if non-aromatic), 4F3, A1IJE, 4JFG-family
}

STANDARD_AA = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
               'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'MSE', 'PHE', 'PRO',
               'SER', 'THR', 'TRP', 'TYR', 'VAL'}


def rotation_matrix_from_vectors(vec1, vec2):
    """
    Find rotation matrix that aligns vec1 to vec2.
    """
    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)
    
    # Check if vectors are already aligned
    if np.allclose(a, b):
        return np.eye(3)
    if np.allclose(a, -b):
        # 180° rotation - need to find an orthogonal axis
        ortho = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        ortho = ortho - np.dot(ortho, a) * a
        ortho = ortho / np.linalg.norm(ortho)
        return Rotation.from_rotvec(np.pi * ortho).as_matrix()
    
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    
    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])
    
    R = np.eye(3) + vx + np.dot(vx, vx) * ((1 - c) / (s ** 2))
    return R


def analyze_structure(cif_file, pdb_id):
    """
    Analyze structure with proper barrel axis alignment.
    """
    try:
        structure = gemmi.read_structure(str(cif_file))
        model = structure[0]
        
        # Find main chain
        best_chain = None
        best_count = 0
        for chain in model:
            count = sum(1 for res in chain if res.name in STANDARD_AA)
            if count > best_count:
                best_count = count
                best_chain = chain
        
        if best_chain is None or best_count < 100:
            return None
        
        # Collect atoms
        ca_atoms = []
        all_atoms = []
        chromophore_atoms = []
        chromophore_type = None
        
        for residue in best_chain:
            for atom in residue:
                if atom.element.name == 'H':
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
            elif (residue.name not in STANDARD_AA and
                  residue.name not in NON_CHROMOPHORE_MODS and
                  50 <= residue.seqid.num <= 90):
                # Scaffold-based fallback: check imidazolinone atoms
                scaffold_count = sum(1 for atom in residue if atom.name in SCAFFOLD_ATOMS)
                if scaffold_count >= 3:
                    for atom in residue:
                        if atom.element.name != 'H':
                            chromophore_atoms.append(np.array([atom.pos.x, atom.pos.y, atom.pos.z]))
                    chromophore_type = residue.name
        
        ca_atoms = np.array(ca_atoms)
        all_atoms = np.array(all_atoms)
        
        if len(ca_atoms) < 100:
            return None
        
        # ============================================================
        # STEP 1: Calculate barrel axis via PCA of Cα atoms
        # ============================================================
        ca_centroid = np.mean(ca_atoms, axis=0)
        ca_centered = ca_atoms - ca_centroid
        
        cov = np.cov(ca_centered.T)
        eigenvalues, eigenvectors = eigh(cov)
        
        # Sort by eigenvalue (largest = barrel axis direction)
        idx = eigenvalues.argsort()[::-1]
        eigenvalues = eigenvalues[idx]
        eigenvectors = eigenvectors[:, idx]
        
        barrel_axis = eigenvectors[:, 0]
        
        # Ensure consistent orientation (arbitrary but consistent)
        if barrel_axis[2] < 0:
            barrel_axis = -barrel_axis
        
        # ============================================================
        # STEP 2: Rotate structure so barrel axis = z-axis
        # ============================================================
        z_axis = np.array([0, 0, 1])
        R = rotation_matrix_from_vectors(barrel_axis, z_axis)
        
        # Apply rotation to all atoms
        all_atoms_rotated = np.dot(all_atoms - ca_centroid, R.T)
        ca_atoms_rotated = np.dot(ca_atoms - ca_centroid, R.T)
        
        if chromophore_atoms:
            chromophore_atoms = np.array(chromophore_atoms)
            chromophore_rotated = np.dot(chromophore_atoms - ca_centroid, R.T)
            chromophore_center = np.mean(chromophore_rotated, axis=0)
        else:
            chromophore_center = None
        
        # ============================================================
        # STEP 3: Center structure (translate so slice is at z=0)
        # ============================================================
        
        # Option: Center on chromophore (if present) or barrel center
        if chromophore_center is not None:
            center_z = chromophore_center[2]
        else:
            # Use barrel center (mean z of all atoms)
            center_z = np.mean(all_atoms_rotated[:, 2])
        
        # Translate so center is at z=0
        all_atoms_final = all_atoms_rotated.copy()
        all_atoms_final[:, 2] -= center_z
        
        # ============================================================
        # STEP 4: Calculate cross-section at z=0 (now perpendicular to barrel)
        # ============================================================
        
        thickness = 2.0
        mask = np.abs(all_atoms_final[:, 2]) <= thickness
        slice_atoms = all_atoms_final[mask]
        
        if len(slice_atoms) < 10:
            return None
        
        # Project to XY plane
        points_2d = slice_atoms[:, :2]
        
        # Convex hull
        try:
            hull = ConvexHull(points_2d)
            area = hull.volume
            perimeter = hull.area
        except:
            return None
        
        # Fit ellipse via covariance
        cov_2d = np.cov(points_2d.T)
        eigvals_2d, eigvecs_2d = eigh(cov_2d)
        eigvals_2d = np.sort(eigvals_2d)[::-1]
        
        # Axes (4 standard deviations encompasses ~95% of points)
        major_axis = 4 * np.sqrt(eigvals_2d[0])
        minor_axis = 4 * np.sqrt(eigvals_2d[1])
        
        # Eccentricity
        if major_axis > 0:
            eccentricity = np.sqrt(1 - (minor_axis / major_axis) ** 2)
        else:
            eccentricity = 0
        
        # Circularity
        if perimeter > 0:
            circularity = 4 * np.pi * area / (perimeter ** 2)
        else:
            circularity = 0
        
        # Calculate barrel length (z-extent after rotation)
        barrel_length = ca_atoms_rotated[:, 2].max() - ca_atoms_rotated[:, 2].min()
        
        # Alignment quality check
        alignment_angle = np.degrees(np.arccos(np.abs(np.dot(barrel_axis, z_axis))))
        
        return {
            'pdb_id': pdb_id,
            'has_chromophore': chromophore_type is not None,
            'chromophore_type': chromophore_type,
            
            # Cross-section metrics
            'n_atoms': len(slice_atoms),
            'convex_area': area,
            'convex_perimeter': perimeter,
            'major_axis': major_axis,
            'minor_axis': minor_axis,
            'eccentricity': eccentricity,
            'circularity': circularity,
            
            # Barrel properties
            'barrel_length': barrel_length,
            'original_axis_angle': alignment_angle,
            
            # Eigenvalues (for shape analysis)
            'eigenvalue_ratio': eigenvalues[0] / eigenvalues[1] if eigenvalues[1] > 0 else 0,
        }
    
    except Exception as e:
        return None


def main():
    print("="*70)
    print("CORRECT CROSS-SECTION ANALYSIS")
    print("Aligned perpendicular to barrel axis")
    print("="*70)
    
    # Find structure directory
    struct_dir = None
    for d in [Path("scop_gfp_structures"),
              Path.home() / "Desktop" / "scop_gfp_structures",
              Path.home() / "Downloads" / "scop_gfp_structures"]:
        if d.exists():
            struct_dir = d
            break
    
    if struct_dir is None:
        print("ERROR: Cannot find scop_gfp_structures directory")
        return
    
    print(f"\nStructure directory: {struct_dir}")
    
    # Process all structures
    cif_files = sorted(struct_dir.glob("*.cif"))
    print(f"Total structures: {len(cif_files)}")
    
    results = []
    failed = []
    
    print("\nProcessing structures...")
    
    for i, cif_file in enumerate(cif_files):
        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(cif_files)}...")
        
        pdb_id = cif_file.stem.upper()
        
        if pdb_id in REMOVE_IDS:
            continue
        
        result = analyze_structure(cif_file, pdb_id)
        
        if result:
            result['expected_chromophore'] = pdb_id not in NO_CHROMOPHORE
            results.append(result)
        else:
            failed.append(pdb_id)
    
    df = pd.DataFrame(results)
    
    print(f"\nSuccessfully analyzed: {len(df)}")
    print(f"Failed: {len(failed)}")
    
    # ================================================================
    # STATISTICS
    # ================================================================
    print(f"\n" + "="*70)
    print("CROSS-SECTION STATISTICS (Properly Aligned)")
    print("="*70)
    
    print(f"\n{'Metric':<20} {'Mean':>12} {'Std':>12} {'Min':>12} {'Max':>12}")
    print("-"*70)
    
    for col in ['convex_area', 'minor_axis', 'major_axis', 'eccentricity', 'circularity', 'n_atoms', 'barrel_length']:
        if col in df.columns:
            print(f"{col:<20} {df[col].mean():>12.2f} {df[col].std():>12.2f} {df[col].min():>12.2f} {df[col].max():>12.2f}")
    
    # ================================================================
    # CHROMOPHORE EFFECT
    # ================================================================
    print(f"\n" + "="*70)
    print("CHROMOPHORE EFFECT (Properly Aligned Cross-Sections)")
    print("="*70)
    
    # Use detected chromophore, not expected
    with_chrom = df[df['has_chromophore'] == True]
    without_chrom = df[df['has_chromophore'] == False]
    
    print(f"\n  WITH chromophore (detected):    n = {len(with_chrom)}")
    print(f"  WITHOUT chromophore (detected): n = {len(without_chrom)}")
    
    # Also check expected
    with_expected = df[df['expected_chromophore'] == True]
    without_expected = df[df['expected_chromophore'] == False]
    print(f"  WITH chromophore (expected):    n = {len(with_expected)}")
    print(f"  WITHOUT chromophore (expected): n = {len(without_expected)}")
    
    print(f"\n{'Metric':<20} {'WITH':>12} {'WITHOUT':>12} {'Δ':>10} {'p-value':>12}")
    print("-"*68)
    
    for col in ['convex_area', 'minor_axis', 'major_axis', 'eccentricity', 'circularity', 'barrel_length']:
        with_vals = with_chrom[col].dropna()
        without_vals = without_chrom[col].dropna()
        
        if len(with_vals) > 5 and len(without_vals) > 5:
            stat, pval = stats.mannwhitneyu(with_vals, without_vals)
            diff = with_vals.mean() - without_vals.mean()
            sig = '***' if pval < 0.001 else ('**' if pval < 0.01 else ('*' if pval < 0.05 else ''))
            print(f"{col:<20} {with_vals.mean():>12.2f} {without_vals.mean():>12.2f} {diff:>+10.2f} {pval:>10.4f} {sig}")
    
    # ================================================================
    # VERIFICATION: Compare to previous methods
    # ================================================================
    print(f"\n" + "="*70)
    print("VERIFICATION: Comparison with Previous Methods")
    print("="*70)
    
    # Load old z=0 data
    old_file = None
    for f in [Path("Results/cross_section_data.csv"), Path("cross_section_data.csv")]:
        if f.exists():
            old_file = f
            break
    
    if old_file:
        df_old = pd.read_csv(old_file)
        df_old['pdb_id'] = df_old['pdb_id'].str.upper()
        
        df_compare = df.merge(df_old[['pdb_id', 'convex_area', 'eccentricity', 'tm_score']], 
                              on='pdb_id', suffixes=('_correct', '_old'))
        
        print(f"\nStructures compared: {len(df_compare)}")
        
        # Correlations
        rho_area, p_area = stats.spearmanr(df_compare['convex_area_correct'], 
                                           df_compare['convex_area_old'])
        rho_ecc, p_ecc = stats.spearmanr(df_compare['eccentricity_correct'], 
                                          df_compare['eccentricity_old'])
        
        print(f"\nCorrelation (Correct vs Old z=0):")
        print(f"  Area: ρ = {rho_area:.3f} (p = {p_area:.4f})")
        print(f"  Eccentricity: ρ = {rho_ecc:.3f} (p = {p_ecc:.4f})")
        
        print(f"\nMean values:")
        print(f"  Area - Correct: {df_compare['convex_area_correct'].mean():.1f}, Old: {df_compare['convex_area_old'].mean():.1f}")
        print(f"  Ecc - Correct: {df_compare['eccentricity_correct'].mean():.3f}, Old: {df_compare['eccentricity_old'].mean():.3f}")
        
        # Does TM-score still correlate with area (it shouldn't now)
        rho_tm, p_tm = stats.spearmanr(df_compare['tm_score'], df_compare['convex_area_correct'])
        print(f"\nTM-score vs Area (Correct method): ρ = {rho_tm:.3f} (p = {p_tm:.4f})")
        
        rho_tm_old, p_tm_old = stats.spearmanr(df_compare['tm_score'], df_compare['convex_area_old'])
        print(f"TM-score vs Area (Old method): ρ = {rho_tm_old:.3f} (p = {p_tm_old:.4f})")
        
        if abs(rho_tm) < abs(rho_tm_old):
            print(f"\n✓ Correct method shows LESS TM-score bias!")
        
    # Save results
    output_file = 'cross_section_correct.csv'
    df.to_csv(output_file, index=False)
    print(f"\nSaved: {output_file}")
    
    # Summary
    print(f"\n" + "="*70)
    print("METHODOLOGY SUMMARY")
    print("="*70)
    print(f"""
WHAT THIS ANALYSIS DOES DIFFERENTLY:

1. BARREL AXIS ALIGNMENT
   - Calculates barrel axis from Cα PCA (not from TM-align)
   - Rotates each structure so barrel axis = z-axis
   - Ensures slice is PERPENDICULAR to barrel

2. CENTERING
   - Centers on chromophore (if present) or barrel center
   - Slice at z=0 is now at the biologically relevant location

3. RESULT
   - True cross-sections, not oblique cuts
   - Comparable across all structures
   - Independent of TM-align reference orientation

EXPECTED DIFFERENCES FROM OLD METHOD:
   - Lower eccentricity (old oblique cuts artificially elongated)
   - More consistent areas (less variation from slice angle)
   - Less TM-score correlation (removes alignment artifact)
""")


if __name__ == "__main__":
    main()
