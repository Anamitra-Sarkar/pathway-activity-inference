# Architecture

## Overview
Context-specific biological pathway activity inference from bulk RNA-seq/microarray expression (samples × genes). The system computes per-sample, per-pathway activity scores using two established single-sample methods, validates them via differential activity testing, and exposes scoring through a FastAPI backend with a React dashboard. Small synthetic fixtures (with honest injected signal) verify correctness in this sandbox; real-run path is documented for Kaggle/Modal with downloaded GMTs and real GEO matrices.

```
Expression Matrix (samples × genes)
      │
      ├─ GMT Parser (MSigDB Hallmark / Reactome) ─┐
      │                                            │
      ├─► ssGSEA (Barbie 2009) ─┐                  │
      │                         ├─► Pathway Scores (samples × pathways)
      └─► z-score (Lee 2008) ───┘                  │
                                                   ├─► Differential (Wilcoxon + BH-FDR)
                                                   └─► Correlation (Spearman/Pearson agreement)
                                                           │
                                ┌──────────────────────────┘
                                ▼
                         FastAPI backend
                         (scoring ungated, curated gated)
                                │
                         React/Vite frontend
                         (heatmap, bar chart, table)
```

## Components

### `data_pipeline/`
- **`gmt_parser.py`** — Parses GMT (pathway_name `\t` description `\t` gene1 `\t` ...). Supports MSigDB Hallmark 50 and Reactome. Handles real-file quirks: UTF-8 BOM (`utf-8-sig`), CRLF/LF, comment lines `#`, blank/whitespace lines, empty gene fields filtered, dedup preserving order, duplicate pathway names raise `ValueError`, empty-file/only-comments raises honest error. Also exposes `parse_gmt_string` with identical strictness. CLI: `python -m data_pipeline.cli --gmt-path <file> --expr <csv>`.
- **`ssgsea.py`** — Real ssGSEA (Barbie et al. Nature 2009; Subramanian et al. PNAS 2005; clarified in Hänzelmann et al. BMC Bioinf 2013).
  - Per sample: rank genes descending by expression, weighting exponent `alpha=0.25` (validated `0 < alpha ≤ 5`), integrated enrichment `ES = sum_i [P_hit(i) - P_miss(i)]` where `P_hit = sum_{j≤i, j∈G} (N - rank_j)^alpha / N_R`, `P_miss = sum_{j≤i, j∉G} 1/(N-|G|)`.
  - Top-ranked gene sets score high positive, bottom low negative, random ~0. Validated in `tests/test_ssgsea.py`. Edge: `k==0` or `k==N` → 0, NaNs dropped per sample, empty matrix raises.
- **`zscore.py`** — Combined z-score (Lee et al. PLoS Comput Biol 2008). Gene-wise z across samples, then per-sample pathway mean (optionally `sum/√k`). Fast baseline to compare against ssGSEA. Handles zero-variance genes (`σ→1`), missing genes → 0.
- **`differential.py`** — Wilcoxon rank-sum (Mann-Whitney U, two-sided, `scipy.stats.mannwhitneyu`) per pathway between two groups, plus Benjamini-Hochberg FDR (1995) across pathways. BH enforces monotonic `q` via reverse cumulative min. Validates labels (2 groups required else 400), missing samples raise honest `ValueError`.
- **`correlation.py`** — Spearman/Pearson per pathway between methods (`scipy.stats.pearsonr/spearmanr`). Honest agreement metric (e.g., injected IFN pathway r>0.5). Validates `method` param (`pearson|spearman|both`), handles constant/n<3 → NaN, no overlapping samples/pathways raises.
- **`cli.py`** — Entry point for real runs with downloaded GMTs and expression CSVs. Writes `ssgsea_scores.csv`, `zscore_scores.csv`, `correlation.csv`, `differential.csv`. Validates `--alpha (0,5]`, file existence, non-empty expression, numeric columns, groups has ≥2 labels, whitespace-trimmed labels.

### `backend/` — FastAPI
- **`main.py`** — Endpoints:
  - `GET /health`, `GET /ready` — honest liveness/readiness; reflect gate truth, never fabricate.
  - `GET /api/v1/pathway-db/status` — curated DB status (gated payload vs honest closed message).
  - `POST /api/v1/score/ssgsea`, `/zscore`, `/both` — **ungated** (deterministic statistical methods, not learned models — documented design choice). Accept `expression: {sample:{gene:value}}`, `gene_sets: {pathway:[genes]}`. Pydantic validation: `expression` non-empty, all values finite numeric, `gene_sets` non-empty, `alpha` in `(0,5]` → 422 on violation; pipeline `ValueError` → 400, unexpected → 500 (never raw 500 with stack trace). `_to_df` coerces numeric and rejects all-NaN.
  - `POST /api/v1/differential` — differential on any score matrix. Validates `scores`/`groups` non-empty, finite; business-logic errors (single group, missing labels) → 400 with clear detail, not 500.
  - `GET /api/v1/curated/status-gated` — example gated artifact endpoint (gate only).
  - `GET /api/v1/curated/analysis` — gated **and** auth-protected (gate + Firebase bearer).
  - CORS enabled, Pydantic schemas for validation with explicit `Field` constraints.
- **`release_gate.py`** — Fail-closed gate: curated endpoints require `MODEL_RELEASE_APPROVED=true` **and** `APPROVED_ARTIFACT_REVISION` non-empty. `gate_status()` honest dict; `require_release_approved` dependency raises 403 with hint when closed. Scoring endpoints deliberately bypass gate (documented).
- **`auth.py`** — Firebase-auth-shaped bearer verification. Reads `FIREBASE_SERVICE_ACCOUNT_JSON` path env; if file+`firebase_admin` present, calls `auth.verify_id_token`. Otherwise fallbacks to JWT structure check (3 parts) for sandbox — never fabricates success without Bearer. Supports `FIREBASE_AUTH_DISABLED=true` for local dev (explicit opt-in, logged, unit-tested with mock). Always returns 401 with honest detail when missing/ malformed.
- **`schemas.py`** — Pydantic models with `Field(gt=0, le=5)` for alpha, `field_validator` for `expression`/`gene_sets`/`scores`/`groups` (non-empty, finite, size caps: 5000 samples, 30000 genes/sample, 1000 pathways).

**Release-gate design rationale:** Core scoring is statistical, not a trained model, so blocking it would break local/offline use. The gate protects *promoted artifacts*: a curated, validated pathway DB revision / precomputed reference bundle that would be served in production only after approval. Docs + health truthfully distinguish scoring-ready vs curated-ready.

### `frontend/` — React + Vite + TypeScript
- Non-boilerplate scientific dashboard: expression paste/upload (long `sample,gene,value` or wide `samples×genes`), gene-set JSON paste (mirrors GMT), ranked pathway table with `q_value` pills, bar chart of `-log10(q)` for top 6 pathways, heatmap of per-sample ssGSEA scores (diverging blue→red, min/max normalized), correlation table (Spearman/Pearson). Uses `fetch` to `POST /api/v1/score/both` + `/differential`.
- Vite proxy `/api`→`:8000`. Build: `npm run build`. Dev: `npm run dev`.
- Design: gradient header, cards, 14px radii, muted text `#64748b`, pill significance (green `q<0.05` else gray), monospace inputs. Accessibility: `skip-link`, `<label>`/`htmlFor` for all textareas/selects, `aria-label`/`aria-describedby`/`aria-invalid`/`aria-busy`/`role="alert"`/`role="status"`/`aria-live`, `role="progressbar"` for -log10 bars, `role="grid"`/`gridcell` + screen-reader table alternative for heatmap, `caption`/`scope="col"` on tables, `focus-visible` outlines, responsive at 900px and 480px (container padding, header size, `overflow-x:auto` tables), `sr-only` utilities.
- Validation: client shows `Invalid JSON` alert with `role="alert"`, loading `aria-busy`, error `aria-live="assertive"`, health status `aria-live="polite"`.

### `tests/` — pytest, synthetic fixtures
- `conftest.make_injected_fixture()` — 4 samples × 13 genes, IFN genes injected high in A (9) vs B (3); background ~5. Honestly synthetic; verifies pipeline recovers known signal (not a biological claim).
- `test_ssgsea` — high vs random vs low ES; shape; injected recovery.
- `test_zscore` — shape, scaling `√k`, injected recovery.
- `test_gmt` — parsing, dedup, malformed, missing file.
- `test_differential` — BH monotonicity/bounds; injected differential top hit; explicit group swapping.
- `test_correlation` — injected agreement `r>0.5`; constant handling.
- `test_api` — health/ready honesty, scoring ungated, gated 403→200 when env set, auth 401/200 paths, `optional_auth` for scoring.
- `test_hardening` — edge-case coverage: GMT comments/BOM/CRLF/whitespace/dedup/empty-name/duplicate variants, backend validation 422 for empty/alpha out-of-range/non-finite/missing fields, differential single-group 400, correlation invalid method, ssgsea alpha bounds, ungated-gate still open, pathway-db honest status.

## Data Flow (real run)
1. Download GMT: `wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt` or Reactome `https://reactome.org/download/current/ReactomePathways.gmt`.
2. Prepare expression CSV (e.g., GEO series like GSE183947, no auth; rows samples, columns genes).
3. Run: `python -m data_pipeline.cli --gmt-path data/h.all.v2023.2.Hs.symbols.gmt --expr data/expression.csv --outdir results/ --groups data/groups.csv` (or `--expr-transpose`).
4. Serve backend: `MODEL_RELEASE_APPROVED=true APPROVED_ARTIFACT_REVISION=$(git rev-parse HEAD) uvicorn backend.main:app --host 0.0.0.0 --port 8000` (or without env: scoring still works, curated blocked).
5. Frontend: `cd frontend && npm install && npm run dev` (proxies to `:8000`).

## Verification (this repo)
- All tests use tiny synthetic matrices (no downloads). Run `pytest -q` — 69 tests: correctness properties, GMT quirks (BOM/comments/CRLF), validation 422/400, gate, auth stub.
- Frontend builds via Vite TS check (`npm run build` typechecks); backend tested via `TestClient`. A11y verified via labels/roles/focus-visible and narrow-width responsive.

## Citations
- Barbie DA et al. Nature 2009 (ssGSEA) — PMID 19648913.
- Subramanian A et al. PNAS 2005 (GSEA) — PMID 16199517.
- Hänzelmann S et al. BMC Bioinf 2013 (GSVA clarifies ssGSEA) — PMID 23323831.
- Lee E et al. PLoS Comput Biol 2008 (combined z-score) — PMID 18989396.
- Benjamini & Hochberg JRSS-B 1995 (FDR).
- Mann & Whitney 1947 / Wilcoxon 1945 (rank-sum).

## Non-Goals / Out-of-Scope
- No heavy compute or downloads in tests/sandbox; real GMT downloads documented with `--gmt-path`.
- No trained ML model registry; gate applies to curated DB revision, not to statistical functions.
- No Firebase Admin SDK required in sandbox; stub tests mock verifier.
