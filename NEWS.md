# RFindsc 1.0.0

- First stable release, the version of the RFind-sc paper. The code and the scores are those of 0.3.9; only the
  version number changed.

# RFindsc 0.3.9

- `score_cells()` and the functions built on it are faster on several cores: the per-cell signatures are sent back from
  the worker processes only with `return_signatures = TRUE` (before, always; for 75,572 cells about 300 MB), and the
  counts are copied into a dense matrix before the workers start only when the copy takes at most 2 GB (before, 8 GB;
  `RFIND_DENSE_MAX_GB` still sets the cap). Scoring 75,572 cells against one panel of 500 genes per side took 46 s on
  6 cores instead of 70 s (R 4.5.2, Apple M3, 16 GB), and the 5.7 GB copy no longer fills the memory of an 8 GB machine.
  Scores are unchanged.
- `running_fisher()` warns (class `rfindsc_duplicated_names`) when a list names a gene more than once, which counts the
  gene twice and inflates the overlap; `score_cells()` warns once when names of the scoring universe differ only in case
  (e.g. PISD and Pisd), since names are matched without regard to case. Scores are unchanged.
- README: the comparison with AUCell and UCell describes UCell's signed signatures (up minus down genes, set to 0 when
  negative), as the RFind-sc paper does.

# RFindsc 0.3.8

- The bundled references (`default_biosets()`: DAM, Maturity, Interferon, G2M; `load_dedonno_biosets()`: Glucocorticoid,
  Seizure, UPR, NF-kB) are rebuilt from the final references of the RFind-sc paper (2026-10-16): counts per million over the
  nuclear genes, genes without informative names, mitochondrial and haemoglobin genes left out and ribosomal protein genes
  kept, the top 500 genes per direction of each published list taken before the genes of a dataset, and the top 200
  up-regulated genes for Seizure without an expression floor. Scores against DAM, Maturity, Glucocorticoid and Seizure
  change accordingly; the gene sets are unchanged.
- `score_cells()` and the functions built on it (`run_pipeline()`, `score_panel_pairs()`): A ranked panel now keeps its top `top_n` genes per side before its genes are restricted to the pair universe, as in the
  paper; before, the lists were restricted first and then cut, so lower-ranked genes of a reference took the place of top
  genes that the dataset did not measure. Gene sets are used whole (before, they were also cut to `top_n` after the
  restriction, which mattered only for sets of more than 500 genes). A panel whose top `top_n` genes per side are all in the
  pair universe gives the same scores as before.

# RFindsc 0.3.7

- New `add_rfind_scores()`: writes the cell x panel scores of `run_pipeline()`, `score_panel_pairs()`,
  `score_genesets()` or `score_cells()` into the metadata of a Seurat object (or a data frame of cells), one column
  per panel (`RFind_<panel>`), matched by cell name, so that they can be plotted with `Seurat::FeaturePlot()`; the
  counterpart of `Seurat::AddModuleScore()`. Scores are not changed. The Python package has the same function for
  `AnnData.obs`.
- `rfindsc_config()`: the default `n_cores` is 1 on Windows, where `parallel::mclapply()` cannot fork (before, scoring with
  the default settings stopped there on a machine with 3 or more cores). Unchanged on other systems.
- `score_cells()`: the cap on the dense copy of the counts made before the workers start (8 GB) can be set with the
  environment variable `RFIND_DENSE_MAX_GB`, as `RFIND_CHUNK_SIZE` sets the chunk size; both are described in
  `?score_cells`. A value of `RFIND_CHUNK_SIZE` that is not a number now gives the default (200), as it does for
  `RFIND_DENSE_MAX_GB` (8); before, scoring stopped with an error from `seq()`. Scores do not depend on either.
- Seurat input needs Seurat >= 5.0.0 (the `layer` argument of `GetAssayData()`): Suggests says so, and scoring a Seurat
  object stops on an older SeuratObject. A split `counts` layer (e.g. after `merge()`) must be joined first
  (`SeuratObject::JoinLayers(x, assay = "RNA")`; without `assay`, only the default assay is joined).
- Suggests no longer lists UCell, homologene and jsonlite, which the package does not use, and asks for msigdbr >= 10.0.0
  (the `collection` / `subcollection` arguments that `msigdbr_to_panels()` and `c2_cgp_auto_pair()` pass; older
  versions name them `category` / `subcategory`). `LICENSE.md` holds the full
  MIT licence text; `data-raw/build_bundled_biosets.R` reads the references from the paths of the RFind-sc paper
  repository (<https://github.com/tmurano/RFind-sc>). Scores unchanged. The bundled `biosets.rds` is the 0.3.5 bundle
  (built 2026-09-18) without DAM2 and Cholesterol, not rebuilt; its `meta$source` still names the earlier location of the
  Hammond file, `Reference/output/references_panel9.rds`, now `data/Hammond2019/output/references/references_panel9.rds` in
  the paper repository.

# RFindsc 0.3.6

- Removed the bundled DAM2 (Sala Frigerio 2019) and Cholesterol (MSigDB Hallmark) references, which the RFind-sc paper
  does not use. `default_biosets()` now returns the 4 Hammond axes (`DAM`, `Maturity`, `Interferon`, `G2M`); their gene
  lists and scores, and those of `load_dedonno_biosets()`, are unchanged.

# RFindsc 0.3.5

- **Pair universe (scoring change for ranked references).** The background of each comparison is now the genes of the
  dataset that the reference could contain, as in the RFind-sc paper (2026-09-18): the scoring universe intersected with the
  panel's new optional field `eligible` (the genes the reference's source study reported with a finite, non-zero change).
  Both the panel's lists and each cell's UP / DOWN lists are restricted to this universe before the top `config$top_n` /
  `config$n_cell` genes are taken. A panel without `eligible` (a membership gene set, which can contain any gene) uses the
  whole scoring universe, as before. Scores of ranked references with `eligible` change (the universe shrinks, so |scores|
  shrink); within one reference all cells share the universe, so the comparison of cells is largely kept.
- The universe size counts distinct lower-case symbols (e.g. `PISD` and `Pisd` count once; before, every row counted). This
  changes scores slightly only when two rows differ only in case.
- `build_panel_from_deg()` records `eligible`: every gene of the table with a finite, non-zero log2FC (before `padj_thresh`
  and `top_n`). Pass the full DEG table (all tested genes) so that it describes what the study measured.
- `score_panel_pairs()` passes a panel's `eligible` through; `to_b2()` keeps it (lower case).
- `score_cells()` returns `universe_size` (the pair universe size per panel).
- Bundled references: `DAM`, `DAM2`, `Maturity`, `Glucocorticoid` and `Seizure` carry `eligible` (from the paper's
  references); the gene groups have none. Gene lists unchanged (`data-raw/build_bundled_biosets.R`).
- New internal helpers: `rank_cell_residuals()` (a cell's full ranking), `prepare_panel()` (pair universe and lists),
  `score_cell_ranking()` (cut inside the pair universe, then the four-term score). `build_cell_signature()` is unchanged.
- Checked against the paper pipeline (`tools/check_rfindsc.R` in the
  [RFind-sc paper repository](https://github.com/tmurano/RFind-sc)): the six Hammond axes and the four De Donno axes
  agree to 1e-13. The shared Python fixture (`tests/testthat/fixtures/py_grouped`) now stores `eligible` and `ranked` so
  that both languages score the same panels.

# RFindsc 0.3.4

- **Gene sets are used whole (scoring fix).** A membership gene set (an MSigDB Hallmark set such as the bundled
  `Interferon`, `G2M`, `Cholesterol`, `UPR` and `NF-kB` groups) has no order, but until 0.3.3 it was passed to the
  Running Fisher kernel as if it were ranked, in the order stored (alphabetical for the bundled groups), so the top-K
  prefixes of the set were scanned too and the score depended on the alphabetical order of gene names. As for a
  biogroup in the original Running Fisher, a set is now used whole: K is fixed at the set size and only the top-G
  cut-offs of the cell's ranked list are scanned. Scores of gene sets change; scores of ranked references (signed DEG
  lists, and UP-only lists ranked by fold change such as `Seizure`) are unchanged.
- `running_fisher()` (and `bid_running_fisher()`) gain `g2_ranked` (default `TRUE`); `g2_ranked = FALSE` gives the
  set form. `rfind_score_cell()` reads `b2$ranked` and applies it to all four terms.
- Panels carry an optional logical `ranked`. `to_b2()` passes it on; a panel without it is a set when its `type` is
  `"gene_group"`, `"msigdb_type_b"` or `"msigdb_cgp_paired"`, and ranked otherwise. `score_genesets()`,
  `msigdbr_to_panels()` and `c2_cgp_auto_pair()` mark their panels `ranked = FALSE`; `score_panel_pairs()` passes a
  panel's `ranked` field through (default ranked).
- Bundled references: every panel has `ranked` (`FALSE` for the five gene groups) and `$meta$ranked`; the type of
  `Seizure` changed from `"gene_group"` to `"ranked_up"` (a ranked UP-only list). Gene lists unchanged
  (`data-raw/build_bundled_biosets.R`).

# RFindsc 0.3.3

- **Per-cell-type (grouped) baseline.** `score_cells()` gains `baseline_groups`: one label per cell (e.g. the cell
  type). Each cell's expected counts then use the mean CPM of the cells of its own group instead of one dataset-wide
  baseline, the rule the RFind-sc paper uses for De Donno et al. (baseline per cell type, both conditions pooled, so
  that cell identity does not dominate the signature). Default `NULL` = one baseline over all cells, unchanged
  behaviour; the scoring universe is still chosen over all cells. Same argument name and semantics as in the Python
  package rfindsc. `run_pipeline()`, `score_genesets()` and `score_panel_pairs()` pass `baseline_groups` and `libsize`
  through to `score_cells()`.
- **Man pages regenerated with roxygen2**, so `man/` matches the R sources again (several pages had been edited or
  written by hand; their correct content was moved into the roxygen comments first). Fixes the unclosed `\item{top_n}`
  in `rfindsc_config.Rd`. `tools::checkRd()` is clean on every page.
- **Bundled references regenerated from the final references of the paper** (`data-raw/build_bundled_biosets.R`).
  The Hammond gene universe now holds 8,960 genes: 90 gene names that had lost their dash in an earlier import are
  restored (`H2.D1` -> `H2-D1`), so the UP sides of `DAM`, `DAM2` and `Maturity` now include MHC genes such as `H2-D1`,
  `H2-K1` and `H2-Q7` (each list is still 500 genes; the DOWN sides and the `Interferon`, `G2M` and `Cholesterol`
  groups are unchanged). The De Donno axes (`Glucocorticoid`, `Seizure`, `UPR`, `NF-kB`) are unchanged.

# RFindsc 0.3.2

- The bundled 50-panel MSigDB Hallmark set and its loader `load_hallmark_biosets()` (shipped since 0.1.0) were
  removed. Gene sets of your own can still be scored with `score_cells()`, `score_genesets()` and
  `score_panel_pairs()`; any MSigDB collection, Hallmark included, can be pulled with
  `msigdbr_to_panels(collection = "H")` (needs the `msigdbr` package).
- The `label` fields of `load_dedonno_biosets()` now match the axis names (`Glucocorticoid`, `Seizure`); in 0.3.1
  they still read `CORT` and `Hyperexcitability`. The gene lists are unchanged.

# RFindsc 0.3.1

- `load_dedonno_biosets()` returns the De Donno axes under the paper's index names: `CORT` is now `Glucocorticoid`
  and `Hyperexcitability` is now `Seizure` (a name should say what the index measures). The gene lists are unchanged.

# RFindsc 0.3.0 (2026-09-15, in development)

* **Running Fisher kernel = the RFind-sc paper definition.** `running_fisher()`
  (new; `bid_running_fisher()` is kept as an alias) returns `-log10` of the
  smallest hypergeometric upper-tail p over every pair of cut-offs, top `G` of
  one list x top `K` of the other. Earlier versions scanned one list against the
  whole other list in both directions and averaged the two p-values. Only the
  rank pairs of shared genes are evaluated (identical to the full G x K scan);
  p on the log scale.
* **Bundled references = the paper's axes.** `default_biosets()` now returns the
  6 Hammond axes (`DAM`, `DAM2`, `Maturity`, `Interferon`, `G2M`,
  `Cholesterol`); the internal `Development` / `Aging` / `LPC` and the
  `Immaturity` axis were removed. New `load_dedonno_biosets()` returns the 4
  De Donno axes (`Glucocorticoid`, `Seizure`, `UPR`, `NF-kB`; named `CORT` and `Hyperexcitability` before 0.3.1). Regenerate with
  `data-raw/build_bundled_biosets.R`.
* `score_cells()` gains `libsize`: pass the total UMI over all genes when
  `counts` has already been restricted to a gene universe (as in the paper).
* `score_cells()` restricts each panel to the scoring universe and then keeps
  the top `config$top_n` genes per side (default 500), as in the paper (before,
  panels were used as given and `top_n` was not applied at scoring).
  `build_panel_from_deg()` therefore keeps all genes by default (`top_n = Inf`),
  so a full-length DEG table gives the paper's reference in any dataset.

# RFindsc 0.2.0 (2026-07-03, in development)

* **Delta method updated to the P2 Pearson residual.** `build_cell_signature()`
  now ranks genes by `(obs - exp)/sqrt(exp + 1)` (observed vs baseline-expected
  counts) instead of the previous log-ratio. The DOWN signature is now every
  gene expressed *below* expectation, ranked by deficit, replacing the old
  dropout-only (CPM == 0) DOWN list — this is robust across cell types and no
  longer biased toward baseline-high housekeeping genes. `score_cells()` passes
  observed counts + expected counts to the signature builder accordingly.
* **Bundled reference swapped to the Hammond 9-axis set** used in the current
  figures: `Development`, `Aging`, `LPC` (Type A), `DAM`, `DAM2`, `Immaturity`
  (Type B), `Interferon`, `G2M`, `Cholesterol` (Type C). `default_biosets()`
  returns these 9 panels. Regenerate with `data-raw/build_bundled_biosets.R`
  from `Reference/output/references_panel9.rds`.

# RFindsc 0.1.0 (2026-04-24, in development)

* Initial skeleton.
* Public API: `run_pipeline()`, `score_cells()`, `bid_running_fisher()`,
  `build_panel_from_deg()`, `default_biosets()`, `msigdbr_to_panels()`,
  `c2_cgp_auto_pair()`, `rfindsc_config()`.
* Bundled: 10-panel microglia reference (Hammond 2019 + external DEG +
  histone).
