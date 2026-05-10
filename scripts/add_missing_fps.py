#!/usr/bin/env python3
"""
Download, analyse, and add the ~30 FP structures missing from the main dataset
due to the SCOP temporal curation gap.

Steps:
  1. Download CIFs from RCSB for the 16 named structures
  2. Query RCSB for additional FP structures from FP-producing organisms
     not already in the 877-structure dataset
  3. Run the same cross-section analysis pipeline
  4. Append to merged_complete_data.csv
"""

import ssl, json, urllib.request, time, sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import ConvexHull
import gemmi

RCSB_CIF = 'https://files.rcsb.org/download/{}.cif'
RCSB_GQL = 'https://data.rcsb.org/graphql'

OUT_DIR  = Path.home() / 'Downloads' / 'deep_analysis'
CIF_DIR  = Path.home() / 'Downloads' / 'missing_fp_structures'
MERGED   = OUT_DIR / 'merged_complete_data.csv'
ADDED    = OUT_DIR / 'added_missing_fps.csv'

CIF_DIR.mkdir(exist_ok=True)

ssl_ctx = ssl.create_default_context()
ssl_ctx.check_hostname = False
ssl_ctx.verify_mode    = ssl.CERT_NONE

STANDARD_AA = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','MSE','PHE','PRO','SER','THR','TRP','TYR','VAL'
}
CHROMOPHORE_RESIDUES = {
    'CRO','CRQ','CR7','CRK','CH7','CH6','CR2','CR8','QYG','GYS',
    'NRQ','5SQ','IIC','HTB','PYG','PYK','KYC','MYK','YCM','ASJ',
    'QMM','SEP','FGL','4IN','4TE','SYG','TYG','SRG','CRF','CFO',
    'CFQ','CFK','CFZ','CR0','CRV','DGR','FGR','FGM','GRB','HCR',
    'I8N','MVR','PGR','RCQ','RCR','RFP','SCR','TFQ','UFG','HHV',
    'KYN','ACD','3TY','HMR','BWV','QKH','M3L','YPR','LPS','B3E',
    'G70','G4P','G3P','0AF','YKY','SHT','KYS','CRS','6FG','1TY',
    'FCO','FPP',
}

SLICE_HALF = 2.0

# ── Exact 31 confirmed missing FPs (provided by Luke) ────────────────────────
ALL_MISSING = {
    '2ZO6': 'Kusabira-Cyan (KCY)',
    '2ZO7': 'Kusabira-Cyan mutant KCY-R1',
    '3BX9': 'mKate (pH 2.0)',
    '3BXA': 'mKate (pH 4.2)',
    '3BXB': 'mKate (pH 7.0)',
    '3BXC': 'mKate (pH 9.0)',
    '3LS3': 'Padron0.9 ON state',
    '3LSA': 'Padron0.9 OFF state',
    '3ST2': 'Dreiklang equilibrium',
    '3ST3': 'Dreiklang OFF state',
    '3ST4': 'Dreiklang ON state',
    '4ZIO': 'mCherry143azF irradiated',
    '5YT1': 'mNeptune684',
    '6S68': 'AausFP2',
    '6U1A': 'FusionRed',
    '6WEM': 'mCrimson 0.9',
    '7DMX': 'PhoCl green form',
    '7DNA': 'PhoCl green+red form',
    '7LQO': 'pnRFP',
    '7LUG': 'pnRFP B30Y mutant',
    '7RHA': 'darkmRuby (pH 5.0)',
    '7RHB': 'darkmRuby (pH 8.0)',
    '7RHC': 'darkmRuby (pH 9.0)',
    '8I4K': 'Azami Red1.0',
    '8PEI': 'SAASoti C21N/V127T',
    '8WGP': 'DsRed-Monomer',
    '8ZBO': 'moxSAASoti F97M',
    '9LD5': 'LSSmOrange cryo pH 8.0',
    '9LD8': 'LSSmOrange RT pH 8.0',
    '9LD9': 'LSSmOrange RT 250ps post-pump',
    '9S0T': 'sfGFP p-(phenylazo)-L-Phe',
}

# Deduplicate against existing dataset
df_existing = pd.read_csv(MERGED)
existing_ids = set(df_existing['pdb_id'].str.upper())

all_candidates = ALL_MISSING

to_process = {k: v for k, v in all_candidates.items() if k not in existing_ids}
print(f'Candidates to process: {len(to_process)}')
print(f'(Already in dataset:   {len(all_candidates) - len(to_process)})')

# ── Download CIF from RCSB ───────────────────────────────────────────────────
def download_cif(pdb_id):
    dest = CIF_DIR / f'{pdb_id.lower()}.cif'
    if dest.exists():
        return dest
    url = RCSB_CIF.format(pdb_id.upper())
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=60, context=ssl_ctx) as r:
            dest.write_bytes(r.read())
        return dest
    except Exception as e:
        print(f'  Download failed {pdb_id}: {e}')
        return None

# ── RCSB metadata fetch ───────────────────────────────────────────────────────
def fetch_rcsb_meta(pdb_id):
    """Get resolution, method, em_max, ex_max from RCSB GraphQL."""
    query = """
    query ($id: String!) {
      entry(entry_id: $id) {
        rcsb_entry_info {
          resolution_combined
          experimental_method
        }
        rcsb_entry_container_identifiers { entry_id }
      }
    }
    """
    try:
        payload = json.dumps({'query': query, 'variables': {'id': pdb_id}}).encode()
        req = urllib.request.Request(
            RCSB_GQL, data=payload,
            headers={'Content-Type': 'application/json', 'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=15, context=ssl_ctx) as r:
            data = json.loads(r.read())
        info = data['data']['entry']['rcsb_entry_info']
        res_list = info.get('resolution_combined') or []
        res = float(res_list[0]) if res_list else None
        method = (info.get('experimental_method') or [''])[0]
        return {'resolution': res, 'method': method}
    except:
        return {'resolution': None, 'method': None}

# ── Cross-section analysis ────────────────────────────────────────────────────
def analyse_cif(cif_path, pdb_id):
    try:
        st = gemmi.read_structure(str(cif_path))
        model = st[0]

        best_chain, best_n = None, 0
        for chain in model:
            n = sum(1 for r in chain if r.name in STANDARD_AA)
            if n > best_n:
                best_n, best_chain = n, chain
        if best_chain is None or best_n < 80:
            return None, 'too few residues'

        ca_list, heavy_list = [], []
        chrom_type = None
        chrom_atoms = []

        for res in best_chain:
            is_std   = res.name in STANDARD_AA
            is_chrom = res.name in CHROMOPHORE_RESIDUES
            if is_std or is_chrom:
                for atom in res:
                    if not atom.is_hydrogen():
                        heavy_list.append([atom.pos.x, atom.pos.y, atom.pos.z])
            if is_std:
                ca = res.find_atom('CA', '\0')
                if ca:
                    ca_list.append([ca.pos.x, ca.pos.y, ca.pos.z])
            if is_chrom and chrom_type is None:
                chrom_type = res.name
                for atom in res:
                    if not atom.is_hydrogen():
                        chrom_atoms.append([atom.pos.x, atom.pos.y, atom.pos.z])

        if len(ca_list) < 80:
            return None, 'too few CA'

        ca_arr    = np.array(ca_list)
        heavy_arr = np.array(heavy_list)

        # PCA barrel axis
        centroid = ca_arr.mean(axis=0)
        cov3d    = np.cov((ca_arr - centroid).T)
        ev3d, evec3d = np.linalg.eigh(cov3d)
        barrel_axis  = evec3d[:, np.argmax(ev3d)]

        zhat = np.array([0.0, 0.0, 1.0])
        v = np.cross(barrel_axis, zhat)
        s = np.linalg.norm(v)
        c = np.dot(barrel_axis, zhat)
        if s < 1e-10:
            R = np.eye(3)
        else:
            vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
            R  = np.eye(3) + vx + vx @ vx * (1 - c) / s**2

        heavy_rot = (heavy_arr - centroid) @ R.T
        ca_rot    = (ca_arr   - centroid) @ R.T
        barrel_length = float(ca_rot[:,2].max() - ca_rot[:,2].min())

        # Chromophore z-level
        if chrom_atoms:
            chrom_arr = np.array(chrom_atoms)
            z_chrom   = float(((chrom_arr - centroid) @ R.T)[:,2].mean())
        else:
            z_chrom = 0.0

        # Slice
        mask  = np.abs(heavy_rot[:,2] - z_chrom) <= SLICE_HALF
        pts2d = heavy_rot[mask,:2]
        if len(pts2d) < 20:
            return None, f'only {len(pts2d)} atoms in slice'

        hull  = ConvexHull(pts2d)
        area  = float(hull.volume)
        perim = float(hull.area)
        circ  = 4 * np.pi * area / perim**2

        cov2d  = np.cov(pts2d.T)
        ev2d,_ = np.linalg.eigh(cov2d)
        ev2d   = np.sort(np.abs(ev2d))[::-1]
        major  = 4 * np.sqrt(ev2d[0])
        minor  = 4 * np.sqrt(ev2d[1])
        ecc    = float(np.sqrt(1 - ev2d[1]/ev2d[0])) if ev2d[0] > 0 else 0.0

        # B-factor ratio (chromophore region vs whole chain)
        all_b, chrom_b = [], []
        for res in best_chain:
            if res.name not in STANDARD_AA: continue
            for atom in res:
                all_b.append(atom.b_iso)
            ca = res.find_atom('CA', '\0')
            if ca:
                pos = (np.array([ca.pos.x, ca.pos.y, ca.pos.z]) - centroid) @ R.T
                if abs(pos[2] - z_chrom) <= SLICE_HALF * 4:
                    for atom in res:
                        chrom_b.append(atom.b_iso)
        b_ratio = np.mean(chrom_b)/np.mean(all_b) if all_b and chrom_b else np.nan

        # Chromophore contacts within 4 Å
        contacts = 0
        if chrom_atoms:
            chrom_np = np.array(chrom_atoms)
            std_atoms = []
            for res in best_chain:
                if res.name not in STANDARD_AA: continue
                for atom in res:
                    if not atom.is_hydrogen():
                        std_atoms.append([atom.pos.x, atom.pos.y, atom.pos.z])
            if std_atoms:
                std_np = np.array(std_atoms)
                for ca in chrom_np:
                    dists = np.linalg.norm(std_np - ca, axis=1)
                    contacts += int((dists <= 4.0).sum())

        return {
            'pdb_id':         pdb_id.upper(),
            'has_chromophore':chrom_type is not None,
            'chromophore_type':chrom_type or '',
            'convex_area':    round(area, 1),
            'convex_perimeter':round(perim, 2),
            'major_axis':     round(major, 2),
            'minor_axis':     round(minor, 2),
            'eccentricity':   round(ecc, 3),
            'circularity':    round(circ, 3),
            'barrel_length':  round(barrel_length, 1),
            'b_factor_ratio': round(b_ratio, 3) if not np.isnan(b_ratio) else None,
            'chrom_contacts': contacts,
            'source':         'added_missing',
            'n_residues':     int(best_n),
        }, None
    except Exception as e:
        return None, str(e)

# ── Main ──────────────────────────────────────────────────────────────────────
print('\n' + '='*65)
print('Adding missing FP structures')
print('='*65)

rows = []
skipped = []

for pdb_id, label in to_process.items():
    print(f'\n{pdb_id} — {label}')

    # Fetch metadata first to filter non-X-ray / low resolution
    meta = fetch_rcsb_meta(pdb_id)
    res  = meta['resolution']
    meth = meta['method']
    print(f'  Resolution={res}  Method={meth}')

    # RCSB returns abbreviated method codes: 'X' = X-RAY DIFFRACTION, 'N' = NMR, 'E' = ELECTRON MICROSCOPY
    if meth and meth.upper() not in ('X', 'X-RAY DIFFRACTION', 'NEUTRON DIFFRACTION'):
        print(f'  SKIP: not X-ray ({meth})')
        skipped.append((pdb_id, label, f'not X-ray: {meth}'))
        continue
    if res and res > 3.5:
        print(f'  SKIP: resolution {res} > 3.5 Å')
        skipped.append((pdb_id, label, f'resolution {res}'))
        continue

    # Download CIF
    cif = download_cif(pdb_id)
    if cif is None:
        skipped.append((pdb_id, label, 'download failed'))
        continue

    # Analyse
    result, err = analyse_cif(cif, pdb_id)
    if result is None:
        print(f'  FAILED: {err}')
        skipped.append((pdb_id, label, err))
        continue

    result['resolution'] = res
    result['match_name'] = label
    print(f'  Area={result["convex_area"]}  Minor={result["minor_axis"]}  '
          f'Major={result["major_axis"]}  Ecc={result["eccentricity"]}  '
          f'Chrom={result["chromophore_type"] or "none"}  BLen={result["barrel_length"]}')
    rows.append(result)
    time.sleep(0.3)

print(f'\n{"="*65}')
print(f'Successfully analysed: {len(rows)}')
print(f'Skipped:               {len(skipped)}')

if skipped:
    print('\nSkipped:')
    for pdb_id, label, reason in skipped:
        print(f'  {pdb_id} ({label}): {reason}')

if rows:
    new_df = pd.DataFrame(rows)
    new_df.to_csv(ADDED, index=False)
    print(f'\nNew rows saved to: {ADDED}')

    # Merge into main dataset
    merged = pd.read_csv(MERGED)
    combined = pd.concat([merged, new_df], ignore_index=True)
    combined.to_csv(MERGED, index=False)
    print(f'Main dataset updated: {len(merged)} → {len(combined)} structures')

    print('\nNew structures summary:')
    print(new_df[['pdb_id','match_name','resolution','chromophore_type',
                  'convex_area','minor_axis','major_axis','eccentricity',
                  'circularity','barrel_length']].to_string(index=False))
