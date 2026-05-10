# GFP Beta-Barrel Geometry Project — Claude Instructions

## What this project is
A systematic computational analysis of fluorescent protein (FP) beta-barrel cross-sectional
geometry across 908 crystal structures from the RCSB PDB. The paper is being submitted to
the Journal of Chemical Information and Modeling (JCIM).

**Key finding:** Red-shifted FPs have narrower, more elliptical barrels. The B-factor ratio
(chromophore/barrel B-factors) is the strongest predictor of quantum yield (ρ = −0.493).
AlphaFold correctly reproduces barrel shape but slightly underestimates cross-sectional
dimensions (~0.5–0.9 Å per axis) due to absence of the mature chromophore.

---

## File locations

### This repo
```
scripts/build_jcim_manuscript.py     ← THE manuscript (run this to regenerate the .docx)
scripts/add_missing_fps.py           ← Downloaded and analysed the 31 missing FP structures
scripts/alphafold_all347.py          ← AlphaFold analysis of all 347 FPbase sequences
scripts/alphafold_wt_comparison.py   ← AlphaFold vs crystal paired comparison
scripts/validate_analysis.py         ← Independent validation of the pipeline (for Maddie)
scripts/gen_supplementary_figs.py    ← Generates all Supporting Information figures
scripts/multiple_testing_corrections.py  ← BH correction across 47 hypothesis tests
scripts/multivariate_analysis.py     ← OLS + Random Forest models
scripts/sensitivity_test.py          ← Slice thickness sensitivity analysis

data/merged_complete_data.csv        ← MASTER dataset: 908 structures × 42 columns
data/Table_S1_all_structures.csv     ← Table S1 for the paper (908 rows, 33 columns)
data/megley_dihedrals.csv            ← τ and φ dihedral angles for 690 structures
data/alphafold_all347.csv            ← AlphaFold results: 327 ok, 17 no model
data/alphafold_vs_crystal_n51.csv    ← Paired AF vs crystal comparison (n=51 proteins)
data/fasta_crystal_matches.csv       ← UniProt→PDB mapping for all 346 FASTA entries
data/FPbase_sequences.fasta          ← Maddie's FPbase UniProt sequence list (347 entries)

figures/pub_figures/                 ← All publication-quality figures (.png)
```

### On Luke's Mac (NOT in this repo — too large for git)
```
~/Downloads/scop_gfp_structures/     ← 908 CIF files (crystal structures)
~/Downloads/missing_fp_structures/   ← 31 additional CIF files
~/Downloads/alphafold_all347_structures/  ← 327 AlphaFold CIF files
~/Downloads/deep_analysis/           ← Output directory for most scripts
~/Downloads/GFP_JCIM_Manuscript.docx ← Current compiled manuscript
```

> **To regenerate CIF files:** Crystal structures download from RCSB automatically when
> you run the analysis scripts. AlphaFold structures download from the EBI API.

---

## Dataset summary
- **908 crystal structures** (877 from SCOP/CATH/Pfam curation + 31 added manually)
- **843** have a validated mature chromophore; **65** do not
- **697** matched to FPbase spectral data (76.8%); **316** have quantum yield values
- **690** have Megley τ/φ dihedral angles computed
- Color class distribution: green (474), red (94), yellow (69), cyan (48), orange (11), blue (1)

## Key columns in merged_complete_data.csv
| Column | Description |
|--------|-------------|
| pdb_id | PDB accession |
| has_chromophore | Boolean — mature chromophore detected |
| chromophore_type | 3-letter residue code (CRO, NRQ, etc.) |
| convex_area | Cross-sectional area (Å²) |
| minor_axis / major_axis | Ellipse axes (Å), from 2D covariance of all slice atoms |
| eccentricity | √(1 − b²/a²) |
| circularity | 4πA/P² |
| barrel_length | End-to-end Cα length (Å) |
| b_factor_ratio | chromophore B-factor / barrel B-factor |
| chrom_contacts | Non-H barrel atoms within 4 Å of chromophore |
| tau_main | Megley τ dihedral (degrees) |
| em_max / ex_max | Emission/excitation maxima from FPbase (nm) |
| lit_qy / lit_ec / lit_brightness | Literature quantum yield, extinction coefficient, brightness |
| resolution | Crystal resolution (Å) |
| color_class | green/red/yellow/cyan/orange/blue |

---

## How the pipeline works
1. **Barrel axis:** PCA of backbone Cα coordinates → PC1 = barrel axis
2. **Rotation:** Rodrigues rotation matrix aligns barrel axis to Z
3. **Slice:** All non-H atoms within ±2.0 Å of the chromophore z-coordinate
4. **Geometry:** Convex hull (area, perimeter) + 2D covariance eigenvalues
   - `major = 4 × √(eigenvalue_max)`, `minor = 4 × √(eigenvalue_min)`
   - This is the **correct formula** — earlier scripts had a bug using hull vertices only
5. **Chromophore level:** Mean z of all chromophore heavy atoms; for AlphaFold, z = 0

> **Important:** The axes use `4 × √(eigenvalue)` of the covariance of **all slice atoms**,
> NOT the hull vertices. An earlier version of alphafold_wt_comparison.py used hull vertices
> with factor 2 — that was a ~10 Å bug that has been fixed.

---

## AlphaFold comparison (n=51 paired)
- 51 proteins from the FPbase FASTA have both an AlphaFold model and a crystal structure
- Matching done via UniProt REST API (https://rest.uniprot.org/uniprotkb/{ID}.json)
- Best-resolution crystal structure selected per protein
- Results: AF faithfully captures shape (eccentricity p=0.21 ns, circularity p=0.66 ns)
  but predicts narrower axes (minor: −0.53 Å p=0.002; major: −0.89 Å p=0.002)
  and longer barrels (+3.7 Å p<0.001) vs crystal structures

---

## Manuscript
The manuscript is built by running:
```bash
python3 scripts/build_jcim_manuscript.py
```
Output: `GFP_JCIM_Manuscript.docx` in the same directory (or `~/Downloads/` if run from there).

The script hard-codes some paths to `~/Downloads/deep_analysis/pub_figures/` for figures.
If running on a new machine, either update those paths or symlink `~/Downloads/deep_analysis/`
to this repo's `data/` and `figures/` directories.

---

## What still needs to be done (as of May 2026)
1. **GitHub repo** — create one and push this directory
2. **Figure for AlphaFold paired comparison** — a paired scatter or Bland-Altman plot
   for the 51-protein AF vs crystal comparison (no figure exists yet)
3. **Table S4** — raw vs BH-corrected p-values for all 47 hypothesis tests
   (referenced in manuscript but CSV not yet generated)
4. **Update manuscript paths** — `build_jcim_manuscript.py` still references
   `~/Downloads/deep_analysis/` for figure files; should be made relative
5. **Update memory files** — `~/.claude/projects/.../memory/project_gfp_barrel.md`
   still says "877 structures" in some places

---

## People
- **Luke Begg** — undergraduate researcher, primary author
- **Maddie** — collaborator, provided the FPbase UniProt sequence list (FASTA)
- **Professor Zimmer** — PI, wrote the advisor checklist

## Python dependencies
```
gemmi, numpy, pandas, scipy, scikit-learn, python-docx, matplotlib, biopython
```
Install with: `pip install gemmi numpy pandas scipy scikit-learn python-docx matplotlib biopython`
