#!/usr/bin/env python3
"""
DETAILED CHROMOPHORE TORSION ANGLE ANALYSIS
============================================
Extract multiple torsion angles that define chromophore geometry:

1. τ (tau) - Main bridge torsion (cis/trans)
2. φ (phi) - Phenol ring twist
3. ψ (psi) - Imidazolinone planarity
4. Additional angles for complete characterization
"""

import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

try:
    import gemmi
except ImportError:
    print("ERROR: pip3 install gemmi")
    exit(1)

# Chromophore residue names
CHROMOPHORE_RESIDUES = {
    'CRO', 'CR2', 'GYS', 'SYG', 'CRQ', 'CRF', 'CRW', 'CRY',
    'NRQ', 'NYG', 'CH6', 'CH7', 'CRG', 'CRU', 'CRV', 'CRS',
    '66A', 'CR0', 'GYC', 'LYG', 'TYG', 'OHD', 'SWG', 'QYG',
    'CR7', 'CR8', 'CR9', 'CRK', 'RC7', 'CFY', 'PIA', 'B2H',
}

REMOVE_IDS = {'9N6A', '7T9X', '6JC5', '9YZQ', '6B0B', '6BBO', '8XPP', '7BWN'}


def calculate_dihedral(p1, p2, p3, p4):
    """Calculate dihedral angle between 4 points in degrees."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    
    # Normalize b2
    b2_norm = b2 / np.linalg.norm(b2)
    
    # Calculate normal vectors
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    
    # Normalize
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    
    if n1_norm < 1e-6 or n2_norm < 1e-6:
        return None
    
    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    
    # Calculate angle
    m1 = np.cross(n1, b2_norm)
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    
    return np.degrees(np.arctan2(y, x))


def get_chromophore_torsions(cif_file):
    """
    Extract all chromophore torsion angles.
    
    GFP chromophore structure (formed from Ser65-Tyr66-Gly67):
    
         O=C2--N2                      (imidazolinone ring)
            |    \
           CA2---C1
            |
           CB2                          (bridge)
            |
           CG2----CD1---CE1            (phenol ring)
                   |      |
                  CD2---CE2
                         |
                        OH (or CZ)
    
    Key torsions:
    - tau (τ): CA2-CB2-CG2-CD1  (bridge, cis/trans)
    - phi (φ): CB2-CG2-CD1-CE1  (phenol twist)
    - psi (ψ): N2-CA2-CB2-CG2   (imidazolinone-bridge)
    - chi (χ): C1-CA2-CB2-CG2   (alternative bridge definition)
    """
    try:
        structure = gemmi.read_structure(str(cif_file))
        model = structure[0]
        
        for chain in model:
            for residue in chain:
                if residue.name in CHROMOPHORE_RESIDUES:
                    # Get all atom positions
                    atoms = {}
                    for atom in residue:
                        atoms[atom.name] = np.array([atom.pos.x, atom.pos.y, atom.pos.z])
                    
                    result = {
                        'chromophore_type': residue.name,
                        'n_atoms': len(atoms),
                    }
                    
                    # Define torsion angle sets to try
                    # Each tuple: (name, atom1, atom2, atom3, atom4, description)
                    torsion_definitions = [
                        # Main bridge torsion (tau) - multiple definitions
                        ('tau_1', 'CA2', 'CB2', 'CG2', 'CD1', 'bridge cis/trans'),
                        ('tau_2', 'CA2', 'CB2', 'CG2', 'CD2', 'bridge cis/trans alt'),
                        ('tau_3', 'C1', 'CA2', 'CB2', 'CG2', 'bridge from imid'),
                        
                        # Phenol ring twist (phi)
                        ('phi_1', 'CB2', 'CG2', 'CD1', 'CE1', 'phenol twist 1'),
                        ('phi_2', 'CB2', 'CG2', 'CD2', 'CE2', 'phenol twist 2'),
                        
                        # Imidazolinone geometry (psi)
                        ('psi_1', 'N2', 'CA2', 'CB2', 'CG2', 'imid-bridge'),
                        ('psi_2', 'C2', 'N2', 'CA2', 'CB2', 'imid planarity'),
                        ('psi_3', 'O2', 'C2', 'N2', 'CA2', 'carbonyl orientation'),
                        
                        # Additional characterization
                        ('chi_1', 'N2', 'C2', 'CA2', 'CB2', 'ring-bridge'),
                        ('chi_2', 'CA2', 'CG2', 'CD1', 'CE1', 'full conjugation'),
                        
                        # Phenol hydroxyl (if present)
                        ('oh_1', 'CD1', 'CE1', 'CZ', 'OH', 'hydroxyl orientation'),
                        ('oh_2', 'CE1', 'CZ', 'OH', 'CG2', 'hydroxyl alt'),
                    ]
                    
                    for name, a1, a2, a3, a4, desc in torsion_definitions:
                        if all(a in atoms for a in [a1, a2, a3, a4]):
                            angle = calculate_dihedral(atoms[a1], atoms[a2], 
                                                       atoms[a3], atoms[a4])
                            if angle is not None:
                                result[name] = angle
                                result[f'{name}_atoms'] = f'{a1}-{a2}-{a3}-{a4}'
                    
                    # Classify configuration based on tau
                    tau = result.get('tau_1') or result.get('tau_2') or result.get('tau_3')
                    if tau is not None:
                        abs_tau = abs(tau)
                        if abs_tau < 30:
                            result['config'] = 'cis'
                        elif abs_tau > 150:
                            result['config'] = 'trans'
                        elif 30 <= abs_tau <= 60 or 120 <= abs_tau <= 150:
                            result['config'] = 'twisted'
                        else:
                            result['config'] = 'intermediate'
                        result['tau_main'] = tau
                    else:
                        result['config'] = 'unknown'
                        result['tau_main'] = None
                    
                    # Calculate planarity (deviation from 0 or 180)
                    if tau is not None:
                        abs_tau = abs(tau)
                        planarity = min(abs_tau, 180 - abs_tau)
                        result['planarity_deviation'] = planarity
                    
                    return result
        
        return None
        
    except Exception as e:
        return None


def main():
    print("="*70)
    print("DETAILED CHROMOPHORE TORSION ANGLE ANALYSIS")
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
    
    # Load spectral data for correlation analysis
    spectral_file = None
    for f in [Path("fpbase_matched_final_v2.csv"), Path("fpbase_matched_final.csv")]:
        if f.exists():
            spectral_file = f
            break
    
    df_spectral = None
    if spectral_file:
        df_spectral = pd.read_csv(spectral_file)
        df_spectral['pdb_id'] = df_spectral['pdb_id'].str.upper()
        print(f"Loaded spectral data: {len(df_spectral)} structures")
    
    # Process structures
    results = []
    
    cif_files = sorted(struct_dir.glob("*.cif"))
    print(f"\nProcessing {len(cif_files)} structures...")
    
    for i, cif_file in enumerate(cif_files):
        if (i + 1) % 100 == 0:
            print(f"  Processed {i+1}/{len(cif_files)}...")
        
        pdb_id = cif_file.stem.upper()
        
        if pdb_id in REMOVE_IDS:
            continue
        
        result = get_chromophore_torsions(cif_file)
        
        if result:
            result['pdb_id'] = pdb_id
            results.append(result)
    
    df = pd.DataFrame(results)
    
    print(f"\nStructures with chromophore torsions: {len(df)}")
    
    # ================================================================
    # TORSION ANGLE STATISTICS
    # ================================================================
    print(f"\n" + "="*70)
    print("TORSION ANGLE STATISTICS")
    print("="*70)
    
    torsion_cols = ['tau_1', 'tau_2', 'tau_3', 'phi_1', 'phi_2', 
                    'psi_1', 'psi_2', 'chi_1']
    
    print(f"\n{'Torsion':<12} {'N':>6} {'Mean':>10} {'Std':>10} {'Min':>10} {'Max':>10}")
    print("-"*60)
    
    for col in torsion_cols:
        if col in df.columns:
            vals = df[col].dropna()
            if len(vals) > 0:
                print(f"{col:<12} {len(vals):>6} {vals.mean():>10.1f} {vals.std():>10.1f} {vals.min():>10.1f} {vals.max():>10.1f}")
    
    # ================================================================
    # CONFIGURATION DISTRIBUTION
    # ================================================================
    print(f"\n" + "="*70)
    print("CHROMOPHORE CONFIGURATION DISTRIBUTION")
    print("="*70)
    
    print(f"\nConfiguration (based on tau):")
    for config in ['cis', 'trans', 'twisted', 'intermediate', 'unknown']:
        n = len(df[df['config'] == config])
        pct = 100 * n / len(df) if len(df) > 0 else 0
        print(f"  {config:<15}: n = {n:>4} ({pct:>5.1f}%)")
    
    # Planarity
    if 'planarity_deviation' in df.columns:
        planar = df['planarity_deviation'].dropna()
        print(f"\nPlanarity deviation from ideal (0° or 180°):")
        print(f"  Mean: {planar.mean():.1f}°")
        print(f"  Std:  {planar.std():.1f}°")
        print(f"  Highly planar (<10°): {len(planar[planar < 10])} structures")
        print(f"  Twisted (>30°): {len(planar[planar > 30])} structures")
    
    # ================================================================
    # CORRELATION WITH EMISSION
    # ================================================================
    print(f"\n" + "="*70)
    print("TORSION ANGLE vs EMISSION WAVELENGTH")
    print("="*70)
    
    if df_spectral is not None:
        df = df.merge(df_spectral[['pdb_id', 'ex_max', 'em_max']], on='pdb_id', how='left')
        
        print(f"\nCorrelations with emission (Spearman):")
        print(f"{'Torsion':<12} {'N':>6} {'ρ':>10} {'p-value':>12}")
        print("-"*42)
        
        for col in torsion_cols + ['tau_main', 'planarity_deviation']:
            if col in df.columns:
                valid = df[[col, 'em_max']].dropna()
                if len(valid) > 20:
                    rho, p = stats.spearmanr(valid[col], valid['em_max'])
                    sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))
                    print(f"{col:<12} {len(valid):>6} {rho:>+10.3f} {p:>10.4f} {sig}")
        
        # Compare emission by configuration
        print(f"\nEmission by configuration:")
        for config in ['cis', 'trans', 'twisted']:
            em = df[df['config'] == config]['em_max'].dropna()
            if len(em) > 3:
                print(f"  {config:<10}: {em.mean():.1f} ± {em.std():.1f} nm (n={len(em)})")
        
        # Statistical test cis vs trans
        cis_em = df[df['config'] == 'cis']['em_max'].dropna()
        trans_em = df[df['config'] == 'trans']['em_max'].dropna()
        if len(cis_em) > 3 and len(trans_em) > 3:
            stat, p = stats.mannwhitneyu(cis_em, trans_em)
            print(f"\n  Cis vs Trans: Δ = {trans_em.mean() - cis_em.mean():.1f} nm, p = {p:.4f} {'*' if p < 0.05 else ''}")
    
    # ================================================================
    # CHROMOPHORE TYPE ANALYSIS
    # ================================================================
    print(f"\n" + "="*70)
    print("CHROMOPHORE TYPE DISTRIBUTION")
    print("="*70)
    
    print(f"\nChromophore types:")
    type_counts = df['chromophore_type'].value_counts()
    for chrom_type, count in type_counts.head(15).items():
        pct = 100 * count / len(df)
        
        # Get mean emission for this type
        em = df[df['chromophore_type'] == chrom_type]['em_max'].dropna()
        em_str = f"{em.mean():.0f} nm" if len(em) > 0 else "N/A"
        
        print(f"  {chrom_type:<6}: n = {count:>4} ({pct:>5.1f}%)  emission: {em_str}")
    
    # ================================================================
    # DETAILED EXAMPLES
    # ================================================================
    print(f"\n" + "="*70)
    print("EXAMPLE STRUCTURES")
    print("="*70)
    
    # Show some cis examples
    cis_examples = df[df['config'] == 'cis'].head(5)
    if len(cis_examples) > 0:
        print(f"\nCis chromophores:")
        for _, row in cis_examples.iterrows():
            em = row.get('em_max', 'N/A')
            em_str = f"{em:.0f}" if pd.notna(em) else "N/A"
            print(f"  {row['pdb_id']}: τ = {row['tau_main']:.1f}°, emission = {em_str} nm")
    
    # Show some trans examples
    trans_examples = df[df['config'] == 'trans'].head(5)
    if len(trans_examples) > 0:
        print(f"\nTrans chromophores:")
        for _, row in trans_examples.iterrows():
            em = row.get('em_max', 'N/A')
            em_str = f"{em:.0f}" if pd.notna(em) else "N/A"
            print(f"  {row['pdb_id']}: τ = {row['tau_main']:.1f}°, emission = {em_str} nm")
    
    # Twisted examples
    twisted_examples = df[df['config'] == 'twisted'].head(5)
    if len(twisted_examples) > 0:
        print(f"\nTwisted chromophores:")
        for _, row in twisted_examples.iterrows():
            em = row.get('em_max', 'N/A')
            em_str = f"{em:.0f}" if pd.notna(em) else "N/A"
            print(f"  {row['pdb_id']}: τ = {row['tau_main']:.1f}°, emission = {em_str} nm")
    
    # Save results
    output_file = 'chromophore_torsions.csv'
    df.to_csv(output_file, index=False)
    print(f"\nSaved: {output_file}")
    
    # Summary
    print(f"\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"""
Torsion angles extracted:
  - tau (τ): Main bridge torsion (cis/trans isomerization)
  - phi (φ): Phenol ring twist  
  - psi (ψ): Imidazolinone-bridge angle
  - chi (χ): Additional ring-bridge angles

Key findings:
  - {len(df[df['config'] == 'trans'])} trans ({100*len(df[df['config'] == 'trans'])/len(df):.0f}%)
  - {len(df[df['config'] == 'cis'])} cis ({100*len(df[df['config'] == 'cis'])/len(df):.0f}%)
  - Cis chromophores emit ~{(df[df['config'] == 'trans']['em_max'].mean() - df[df['config'] == 'cis']['em_max'].mean()):.0f} nm bluer than trans
""")


if __name__ == "__main__":
    main()
