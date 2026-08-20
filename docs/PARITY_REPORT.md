# PARITY_REPORT.md

Comparison of the refactored pipeline's output against the reference
figures. Every number below came from an actual run of the refactored
scripts on the real data.

## 1. Survival (`08_survival_km.R` vs. `KaplanMeier_TCGA_Nunes.pdf`)

Group sizes — exact match: TCGA all tumors 139 high / 139 low (n=554);
TCGA WNT-wt 35/35 (n=139); Nunes all tumors 266/266 (n=1063); Nunes WNT-wt
75/75 (n=299).

Log-rank p-values — exact match, all 8 panels (raw follow-up and
5yr/10yr-truncated): TCGA all p=5e-04 (p<0.0001 truncated); TCGA WNT-wt
p=0.028 (both); Nunes all p=0.098 (both); Nunes WNT-wt p=0.15 (both).
Visual overlay confirms identical curve shapes, censor marks, axis ranges.

## 2. FBXW7×WNT mutation boxplot (`09_mutation_expression_boxplots.R` vs. `Nunes.FBXW7.pdf`)

Group sizes exact match: FBXW7-wt/WNT-mut n=623, FBXW7-mut/WNT-mut n=79,
FBXW7-wt/WNT-wt n=129, FBXW7-mut/WNT-wt n=9. P-values exact match: HES1
0.0027/0.12, HES5 0.013/0.53, HES6 0.86/0.39. Both point-style pages
reproduced.

## 3. scRNA-seq QC (`02_preprocess_scrna.py`)

85,571 cells x 26,628 genes — exact match to the manuscript. (Bug caught
and fixed en route: the metaprogram-removal regex was accidentally
case-insensitive and stripped ~1,460 real genes, e.g. PTEN, PIK3CA, in
addition to the intended `pEpi*` pseudo-genes. Fixed to the case-sensitive
`^p[A-Z]` pattern the original notebooks used.)

## 4. CXCR4 stratification (`05_cxcr4_stratification.py`)

Exact match: CXCR4-high 15 patients / 12,903 cells, CXCR4-low 12 / 18,342,
CXCR4-intermediate 35 / 54,326. `Barplots_CXCR4_pos_neg.png` visually
pixel-identical to the reference. Signature-scoring "genes found"
numerators match the original notebooks' cached output exactly
(Notch_Sander 14/14, Fetal 265, RSC 214, CBC 325, wnt_msigdb 124/131,
Fetal_curated 5, WNT_curated 5/5); denominators differ by 1-2 because
this release's signature TSVs exclude an Excel title-row artifact and one
corrupted gene symbol present in the original spreadsheet reads — since
those excluded entries never matched a real gene symbol either way, the
actual scores are numerically identical.

## 5. UMAP figures (`06_dimreduction_umap.py`)

`UMAP_CXCR4_high_low_patients.pdf` — exact match, including the CXCR4
colorbar's `vmax 90th pct: 1.51` value to two decimal places, and
identical cluster shapes/positions/per-patient colors.

`UMAP_3Doners.pdf` — two bugs caught and fixed via direct visual diff
against the reference: (1) RSC was colored green instead of purple
(wrong signature/color grouping); (2) point size was `s=1` (correct for
the 31,245-cell combined plot) instead of the `s=25` the actual source
notebook (`UMAP_3Doners.ipynb`) uses for its much smaller per-donor cell
counts, which made every panel look sparse compared to the reference.
After both fixes: CXCR4 colorbar `vmax 90th pct` matches exactly per donor
(C160: 1.25, C172: 1.05, C168: 1.25), point density matches the
reference, cluster shapes are visually consistent. Cell counts per donor
exact: C160=2013, C172=341, C168=693.

`UMAP_CXCR4_high_low.pdf` (an 11-panel gene/signature grid) was removed
from the pipeline — it was never one of the five requested figures.

## 6. Palantir/CellRank trajectory (`07_trajectory_palantir.py` vs. `Palantir_3Donors_Grid.pdf`)

Root cell (`C168_T_1_1_0_c1_v3_id-GACACGCCACCACTGG`) and terminal cell
(`C168_T_1_1_0_c1_v3_id-TCAGTCCCAAAGCTAA`) are byte-for-byte identical to
the original notebook's cached output, as is the pseudotime range
(0.0–1.0) and cell count (n=693). The "Palantir trajectory" panel is
visually pixel-identical to the reference.

Three bugs caught and fixed:
- Windows multiprocessing crash — scikit-learn's KMeans (terminal-cell
  selection) spawns joblib/loky worker processes, which on Windows
  requires the script's entry point guarded by `if __name__ == "__main__":`.
- Wrong gene panel — an earlier draft used the gene list from the older,
  separate `Palantir_CellRank_3Donors.ipynb` instead of this figure's
  actual source notebook's own `GENE_PANEL` (CXCR4, LGR5, AXIN2, NKD1,
  ANXA3, ANXA1, TACSTD2).
- Wrong panel structure — an earlier draft simplified the original 4-row
  grid (trajectory / gene-trend heatmap / branch expression / CellRank
  projection) into a 3-panel line-plot version. Rebuilt to reproduce the
  original notebook's exact render functions (`plot_gene_trend_heatmaps`,
  the min-max-normalized 3-gene branch-expression plot, and the CellRank
  streamline projection with the correct Spectral_r colormap), collapsed
  to 1 column since only one donor runs.

After all three fixes, all 4 panels match the reference: identical
heatmap colors/scale, identical branch-expression curve shapes/colors/
timing, and the correct rainbow CellRank colormap (previously rendered as
a plain pink/white sequential palette due to the simplified panel).

## Summary

Every group size, p-value, and gene-scoring numerator checked came out an
exact match. Bugs caught and fixed during this build, all via direct
comparison against the real reference figures or manuscript-stated
numbers — not hypothetical:
1. Case-insensitive metaprogram-removal regex (stripped real genes)
2. `collapse_position_to_gene()` losing matrix dimensions on 1-column input
3. `classify_fbxw7_wnt_groups()` not preserving sample names
4. `wnt_msigdb` key mismatch between config and plotting code
5. Missing PCA/neighbors distinction (20 computed vs. 10 used per-patient)
6. RSC miscolored green instead of purple in the 3-donor grid
7. Wrong point size in the 3-donor grid (s=1 instead of s=25)
8. Windows multiprocessing crash in Palantir's terminal-cell KMeans
9. Wrong gene panel in the trajectory figure
10. Oversimplified trajectory figure (3 panels instead of the source's 4)
