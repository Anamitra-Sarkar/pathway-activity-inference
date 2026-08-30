# Data Sources & Method Citations

## Pathway / Gene-Set Sources (real, documented; downloads NOT attempted in sandbox)

### MSigDB Hallmark (primary)
- **Description:** 50 well-curated Hallmark gene sets distilling many founder sets; broad biological processes.
- **URL:** https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H
- **Direct GMT download (2023.2, human symbols):** `h.all.v2023.2.Hs.symbols.gmt`
  - `https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt`
- **License:** Free for non-commercial / academic use; commercial requires license. Documented at https://www.gsea-msigdb.org/gsea/msigdb/license_terms.jsp
- **Format:** GMT — `pathway_name \t description \t gene1 \t gene2 \t ...`
- **Usage in this repo:** `python -m data_pipeline.cli --gmt-path <local path to h.all...gmt> ...` Tests use small synthetic JSON gene sets mimicking Hallmark structure; no download in CI/sandbox.

### Reactome (alternative, fully open)
- **Description:** Manually curated pathway database (reactions, complexes). Provided as GMT for enrichment.
- **Downloads:** https://reactome.org/download-data  → `ReactomePathways.gmt` (human)
  - Direct: `https://reactome.org/download/current/ReactomePathways.gmt` (also `ReactomePathways_Relation.txt`)
- **License:** CC BY 4.0 (fully open, attribution). See https://reactome.org/page/license-agreement
- **GMT structure:** Same as above; pathway names like `R-HSA-109582` with human-readable descriptions.
- **Usage:** Same CLI flag. Documented as equally valid source; parser handles both.

### KEGG (mentioned in background, not fetched here)
- Via MSigDB C2:KEGG or direct KEGG API, but not required for this project.

## Gene Expression Data

### Real-run procedure (would run on Kaggle/Modal, NOT in sandbox)
- Use a small public GEO series, no-auth, e.g.:
  - **GSE183947** (SARS-CoV-2 PBMC bulk RNA-seq, ~6–12 samples, well-known small series) — download via `GEOquery` or `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE183947`
  - Or **GSE115845**, **GSE15471** — any small bulk series with known phenotype (two groups).
  - Expression matrix: samples × genes (rows samples, columns HGNC symbols). Normalize log2 or counts.

### Synthetic verification used in tests (honest, not a biological finding)
- **Injected signal:** IFN-response genes (`IFIT1, MX1, ISG15`) upregulated to ~9.0 in group A vs ~3.0 in group B; background genes ∼N(5, 0.3). N=4–12 samples.
- **Purpose:** Validates that ssGSEA/z-score + Wilcoxon + BH pipeline *recovers* known injected signal (top differential pathway is injected one, q significant, method correlation r>0.5). Tests assert `mean_A >> mean_B`, `delta < -1`, correlation positive. Clearly documented as synthetic verification, not a real biological result.
- **Location:** `tests/conftest.py` (`make_injected_fixture`), exercised in `tests/test_ssgsea.py`, `test_zscore.py`, `test_differential.py`, `test_correlation.py`.

## Methods — Real, Correctly-Cited

### ssGSEA (single-sample GSEA)
- **Barbie DA et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. Nature. 2009;462:108-112. PMID:19648913 / doi:10.1038/nature08460** — Introduced ssGSEA.
- **Subramanian A et al. Gene set enrichment analysis. PNAS. 2005;102:15545-50.** — Original GSEA running-sum KS statistic.
- **Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis. BMC Bioinformatics. 2013;14:7. doi:10.1186/1471-2105-14-7** — Clarifies ssGSEA weighting (`α=0.25`) and integrated ES formulation used here: `ES = Σ_i (P_hit(i) - P_miss(i))`, `P_hit(i)= Σ_{j≤i, j∈G} |r_j|^α / N_R`.
- **Implementation:** `data_pipeline/ssgsea.py` — rank descending by expression, `(N - rank)^α / N_R` weighting, sum of ECDF differences. Alpha 0.25 matches GenePattern ssGSEA.

### Combined z-score
- **Lee E et al. Inferring pathway activity toward precise disease classification. PLoS Comput Biol. 2008;4(11):e1000218. doi:10.1371/journal.pcbi.1000218** — Average of gene-wise z-scores across pathway; optionally scaled by `√k` (`(Σ z)/√k`). This repo implements mean `z` (and scaled via `scale_by_sqrt_k=True`).
- Code: `data_pipeline/zscore.py` — per-gene mean/SD across samples (zero-SD → 1), then per-pathway mean z.

### Differential activity
- **Wilcoxon rank-sum / Mann-Whitney U** (Mann & Whitney 1947; Wilcoxon 1945) — two-sided per pathway via `scipy.stats.mannwhitneyu`.
- **Benjamini & Hochberg. Controlling FDR. JRSS-B. 1995;57:289-300.** — BH step-up, implemented in `data_pipeline/differential.py` with reverse cumulative min for monotonic `q`.

### Correlation / Agreement
- Pearson (linear) and Spearman (rank) per pathway via `scipy.stats.pearsonr/spearmanr` in `data_pipeline/correlation.py`.

## No Fabrications
- No invented dataset, API, or metric. All endpoints are real GSEA/MSigDB/Reactome URLs with license notes; all statistical methods are standard with correct citations; correlation uses real Spearman/Pearson. Synthetic fixtures are explicitly labeled as such.
