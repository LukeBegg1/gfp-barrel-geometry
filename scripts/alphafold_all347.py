#!/usr/bin/env python3
"""
AlphaFold cross-section analysis for all 347 FP sequences in the FASTA.

Downloads AlphaFold structures, runs PCA-based barrel cross-section at z=0
(barrel centroid), using the CORRECT axis method matching the crystal pipeline:
  - covariance of ALL heavy atoms in the ±2 Å slice
  - major = 4 * sqrt(eigenvalue_max)
  - minor = 4 * sqrt(eigenvalue_min)

Outputs: alphafold_all347.csv (one row per sequence, with spectral data from FPbase)
Saves progress incrementally — safe to interrupt and restart.
"""

import ssl, json, urllib.request, time, sys
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.spatial import ConvexHull
import gemmi

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTA    = Path.home() / 'Downloads' / 'UniProt_FPbase_List_4.16.26.fasta'
AF_DIR   = Path.home() / 'Downloads' / 'alphafold_all347_structures'
OUT_DIR  = Path.home() / 'Downloads' / 'deep_analysis'
OUT_CSV  = OUT_DIR / 'alphafold_all347.csv'

AF_DIR.mkdir(parents=True, exist_ok=True)

ssl_ctx = ssl.create_default_context()
ssl_ctx.check_hostname = False
ssl_ctx.verify_mode    = ssl.CERT_NONE

STANDARD_AA = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
    'HIS','ILE','LEU','LYS','MET','MSE','PHE','PRO',
    'SER','THR','TRP','TYR','VAL'
}

# UniProt IDs for the 13 FPbase-only headers (from FPbase API lookup)
FPBASE_TO_UNIPROT = {
    'dis3GFP':    'Q8T6T8',
    'mcavGFP':    'Q8T5F2',
    'mc3':        'Q7Z0W7',
    'Fpmcavgr7.7':'Q8MU47',
    'plamGFP':    'B6CTZ5',
    'Fpaagar':    'Q8MMA1',
    'cFP484':     'Q9U6Y3',
    'h2-3':       '',         # no UniProt ID — skip
    'Fpag_frag':  'Q8MMA2',
    'LanYFP':     'B1PNC0',   # = blFP-Y3, PDB 5LTQ
    'CpYGFP':     'Q2MHN7',   # PDB 2DD7
    'laesGFP':    'Q6WV11',
    'pdae1GFP':   'Q6WV08',
}

# ── Parse FASTA → list of (uniprot_id, display_label) ────────────────────────
def parse_fasta(fasta_path):
    entries = []
    with open(fasta_path) as f:
        for line in f:
            if not line.startswith('>'):
                continue
            header = line.strip()[1:]
            if header.startswith('sp|') or header.startswith('tr|'):
                parts = header.split('|')
                uid   = parts[1]
                label = parts[2].split(' ')[0] if len(parts) > 2 else uid
            else:
                uid   = FPBASE_TO_UNIPROT.get(header, '')
                label = header
            if uid:
                entries.append((uid, label))
            else:
                print(f'  Skipping (no UniProt): {header}')
    return entries

# ── Fetch FPbase spectral data for all proteins ───────────────────────────────
def fetch_fpbase_data():
    """Return dict: uniprot_id → {em_max, ex_max, qy, name, has_pdb}"""
    url = 'https://www.fpbase.org/api/proteins/?format=json&limit=2000'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    try:
        with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as r:
            data = json.loads(r.read())
    except Exception as e:
        print(f'FPbase API error: {e}')
        return {}

    lookup = {}
    for p in data:
        uid = p.get('uniprot', '') or ''
        if not uid:
            continue
        name    = p.get('name', '')
        pdbs    = p.get('pdb', []) or []
        states  = p.get('states', []) or []
        # Collect spectral data from all states
        em_vals, ex_vals, qy_vals = [], [], []
        for st in states:
            em = st.get('em_max')
            ex = st.get('ex_max')
            qy = st.get('qy')
            if em: em_vals.append(em)
            if ex: ex_vals.append(ex)
            if qy is not None: qy_vals.append(qy)
        lookup[uid] = {
            'fpbase_name': name,
            'em_max':  min(em_vals) if em_vals else None,   # primary (shortest) em
            'ex_max':  min(ex_vals) if ex_vals else None,
            'qy':      qy_vals[0]  if qy_vals else None,
            'has_pdb': len(pdbs) > 0,
            'pdb_ids': ','.join(pdbs),
        }
    return lookup

# ── AlphaFold download ────────────────────────────────────────────────────────
def download_af(uniprot_id):
    existing = sorted(AF_DIR.glob(f'AF-{uniprot_id}-F1-model_v*.cif'))
    if existing:
        return existing[-1]
    api_url = f'https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}'
    try:
        req = urllib.request.Request(
            api_url, headers={'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json'})
        with urllib.request.urlopen(req, timeout=30, context=ssl_ctx) as r:
            data = json.loads(r.read())
        entry   = data[0] if isinstance(data, list) else data
        cif_url = entry.get('cifUrl')
        version = entry.get('latestVersion', '?')
        if not cif_url:
            return None
    except Exception as e:
        return None
    dest = AF_DIR / f'AF-{uniprot_id}-F1-model_v{version}.cif'
    try:
        req2 = urllib.request.Request(cif_url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req2, timeout=120, context=ssl_ctx) as r:
            dest.write_bytes(r.read())
        return dest
    except Exception:
        return None

# ── Cross-section analysis ────────────────────────────────────────────────────
SLICE_HALF = 2.0
MIN_ATOMS  = 20

def analyse(cif_path):
    """
    Returns dict of cross-section metrics or None.
    Uses CORRECTED method: cov of ALL slice atoms, axes = 4*sqrt(eigenvalue).
    """
    try:
        s = gemmi.read_structure(str(cif_path))
        m = s[0]
        # Longest standard-AA chain
        best_chain, best_n = None, 0
        for chain in m:
            n = sum(1 for r in chain if r.name in STANDARD_AA)
            if n > best_n:
                best_n, best_chain = n, chain
        if best_chain is None or best_n < 80:
            return None

        ca_list, all_coords, all_b = [], [], []
        for res in best_chain:
            if res.name not in STANDARD_AA:
                continue
            ca = res.find_atom('CA', '\0')
            if ca:
                ca_list.append([ca.pos.x, ca.pos.y, ca.pos.z])
            for atom in res:
                if not atom.is_hydrogen():
                    all_coords.append([atom.pos.x, atom.pos.y, atom.pos.z])
                    all_b.append(atom.b_iso)

        if len(ca_list) < 80:
            return None

        ca_coords  = np.array(ca_list)
        all_coords = np.array(all_coords)

        # PCA barrel axis
        centroid = ca_coords.mean(axis=0)
        centred  = ca_coords - centroid
        cov      = np.cov(centred.T)
        vals, vecs = np.linalg.eigh(cov)
        barrel_axis = vecs[:, np.argmax(vals)]
        z = np.array([0, 0, 1.0])
        v = np.cross(barrel_axis, z)
        s_norm = np.linalg.norm(v)
        c = np.dot(barrel_axis, z)
        if s_norm < 1e-10:
            R = np.eye(3)
        else:
            vx = np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])
            R  = np.eye(3) + vx + vx @ vx * (1 - c) / (s_norm ** 2)

        all_rot = (all_coords - centroid) @ R.T
        ca_rot  = (ca_coords  - centroid) @ R.T

        barrel_length = float(ca_rot[:, 2].max() - ca_rot[:, 2].min())

        # Slice at z=0 (barrel centroid) ± SLICE_HALF
        mask = np.abs(all_rot[:, 2]) <= SLICE_HALF
        pts2d = all_rot[mask, :2]
        if len(pts2d) < MIN_ATOMS:
            return None

        # Convex hull → area, perimeter
        try:
            hull = ConvexHull(pts2d)
        except Exception:
            return None
        area  = hull.volume
        perim = hull.area

        # Axes — CORRECT method (all slice atoms, 4*sqrt)
        cov2d = np.cov(pts2d.T)
        evals, _ = np.linalg.eigh(cov2d)
        evals = np.sort(np.abs(evals))[::-1]
        major = 4 * np.sqrt(evals[0])
        minor = 4 * np.sqrt(evals[1])
        ecc   = np.sqrt(1 - evals[1] / evals[0]) if evals[0] > 0 else 0.0
        circ  = 4 * np.pi * area / perim**2 if perim > 0 else 0.0

        # pLDDT (stored as B-factor in AlphaFold CIFs)
        mean_plddt = float(np.mean(all_b)) if all_b else np.nan

        return {
            'area':           round(area,  1),
            'perimeter':      round(perim, 2),
            'major_axis':     round(major, 2),
            'minor_axis':     round(minor, 2),
            'eccentricity':   round(ecc,   3),
            'circularity':    round(circ,  3),
            'barrel_length':  round(barrel_length, 1),
            'mean_plddt':     round(mean_plddt, 1),
            'n_slice_atoms':  int(mask.sum()),
            'n_residues':     int(best_n),
        }
    except Exception as e:
        return None

# ── Main ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print('=' * 70)
    print(f'AlphaFold cross-section analysis — all 347 FP sequences')
    print('=' * 70)

    # Load existing results (resume if interrupted)
    done = set()
    if OUT_CSV.exists():
        existing_df = pd.read_csv(OUT_CSV)
        done = set(existing_df['uniprot_id'].tolist())
        print(f'\nResuming: {len(done)} already done.')

    # Parse FASTA
    entries = parse_fasta(FASTA)
    print(f'Sequences to process: {len(entries)} (skipping {346 - len(entries) + 1} with no UniProt)')

    # Deduplicate (some UniProt IDs appear more than once)
    seen_uid = set()
    unique_entries = []
    for uid, label in entries:
        if uid not in seen_uid:
            seen_uid.add(uid)
            unique_entries.append((uid, label))
    print(f'Unique UniProt IDs: {len(unique_entries)}')

    # Fetch FPbase spectral data
    print('\nFetching FPbase spectral data...')
    fpbase = fetch_fpbase_data()
    print(f'  FPbase entries loaded: {len(fpbase)}')

    # Download + analyse
    todo = [(uid, lbl) for uid, lbl in unique_entries if uid not in done]
    print(f'\nTo download/analyse: {len(todo)}\n')

    rows = []
    n_ok = 0; n_fail_dl = 0; n_fail_an = 0

    for i, (uid, label) in enumerate(todo):
        sys.stdout.write(f'\r  [{i+1:3d}/{len(todo)}] {uid:12s} {label[:30]:30s}')
        sys.stdout.flush()

        # Download
        cif_path = download_af(uid)
        if cif_path is None:
            sys.stdout.write(f'  NO AF MODEL\n')
            row = {'uniprot_id': uid, 'label': label, 'status': 'no_af_model'}
            rows.append(row)
            n_fail_dl += 1
            time.sleep(0.2)
            continue

        # Analyse
        res = analyse(cif_path)
        if res is None:
            sys.stdout.write(f'  ANALYSIS FAILED\n')
            row = {'uniprot_id': uid, 'label': label, 'status': 'analysis_failed',
                   'cif_file': cif_path.name}
            rows.append(row)
            n_fail_an += 1
            continue

        # Spectral data
        fp = fpbase.get(uid, {})
        row = {
            'uniprot_id':    uid,
            'label':         label,
            'fpbase_name':   fp.get('fpbase_name', ''),
            'status':        'ok',
            'has_pdb':       fp.get('has_pdb', False),
            'pdb_ids':       fp.get('pdb_ids', ''),
            'em_max':        fp.get('em_max'),
            'ex_max':        fp.get('ex_max'),
            'qy':            fp.get('qy'),
            'cif_file':      cif_path.name,
            **res,
        }
        rows.append(row)
        n_ok += 1
        sys.stdout.write(f'  Area={res["area"]:.0f} Minor={res["minor_axis"]:.1f} '
                         f'Major={res["major_axis"]:.1f} Ecc={res["eccentricity"]:.3f} '
                         f'pLDDT={res["mean_plddt"]:.1f}\n')

        # Append to CSV incrementally
        new_df = pd.DataFrame(rows)
        if OUT_CSV.exists():
            old_df = pd.read_csv(OUT_CSV)
            combined = pd.concat([old_df, new_df], ignore_index=True)
        else:
            combined = new_df
        combined.to_csv(OUT_CSV, index=False)
        rows = []   # clear after saving

        time.sleep(0.15)   # polite rate-limit

    print(f'\n\nDone.')
    print(f'  OK:             {n_ok}')
    print(f'  No AF model:    {n_fail_dl}')
    print(f'  Analysis fail:  {n_fail_an}')
    print(f'\nResults saved to: {OUT_CSV}')

    # Summary statistics
    final_df = pd.read_csv(OUT_CSV)
    ok = final_df[final_df['status'] == 'ok']
    print(f'\nGeometry summary (n={len(ok)}):')
    for col in ['area','minor_axis','major_axis','eccentricity','circularity','barrel_length']:
        if col in ok.columns:
            print(f'  {col:16s}: mean={ok[col].mean():.2f}  SD={ok[col].std():.2f}  '
                  f'min={ok[col].min():.2f}  max={ok[col].max():.2f}')
