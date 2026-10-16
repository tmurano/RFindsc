# RFindsc

**Direction-aware panel scoring for single-cell RNA-seq.**

RFind-sc scores each cell against *case-vs-control* gene panels with the
Running Fisher test: the cell's ranked gene list and a reference list are
compared by the smallest hypergeometric p-value over every pair of cut-offs
(top G genes of one list × top K of the other; an unordered gene set such as an
MSigDB Hallmark set is used whole, and only the cell's cut-offs are scanned), and
the four up/down pairs are combined into one **signed** score — positive for case-aligned cells, negative
for control-aligned — which directionless tools like AUCell and UCell cannot
express by design.

Each cell's own UP/DOWN signature is ranked by the **P2 Pearson residual**
`(obs - exp)/sqrt(exp + 1)` of observed vs baseline-expected counts — a
"delta" whose DOWN side is every gene expressed *below* expectation (not just
dropouts).

Version 0.3.8 matches the RFind-sc paper (Murano et al., in preparation):
G x K Running Fisher kernel (gene groups used whole, 0.3.4), the pair universe
of each comparison (the dataset's genes that the reference could contain, 0.3.5),
the paper's reference axes and the paper's per-cell-type baseline option.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("tmurano/RFindsc@v0.3.8")   # the version of the paper; without @v0.3.8, the current version
```

Required dependencies: `Matrix`, `parallel` (base R).
Optional (used by specific functions):
`Seurat` (>= 5.0.0; for `run_pipeline()`, `score_genesets()` and `score_panel_pairs()` on Seurat objects),
`msigdbr` (>= 10.0.0, which has the `collection` / `subcollection` arguments; for `msigdbr_to_panels()` /
`c2_cgp_auto_pair()`).

```r
# Optional extras
install.packages(c("Seurat", "msigdbr"))
```

Scoring runs in parallel with `parallel::mclapply()`. Windows cannot fork, so there the default is one core
(`n_cores = 1`, see [Configuration](#configuration)); scores are the same.

**Python:** [rfindsc-py](https://github.com/tmurano/rfindsc-py) is the Python implementation (same algorithm, same
bundled references, numerically identical scores; AnnData input):

```bash
pip install git+https://github.com/tmurano/rfindsc-py@v0.3.8
```

---

## Quick start

```r
library(RFindsc)
library(Seurat)

mg  <- readRDS("my_microglia_seurat.rds")

# Default: the 4 Hammond axes of the paper
# (DAM / Maturity signed references + Interferon / G2M gene groups)
out <- run_pipeline(mg, biosets = default_biosets())

head(out$scores)   # cell x panel matrix: DAM, Maturity, Interferon, G2M

# Write the scores into the metadata (columns RFind_DAM, RFind_Maturity, ...) and plot them
mg <- add_rfind_scores(mg, out)
FeaturePlot(mg, c("RFind_DAM", "RFind_Maturity"))
```

A signed score matrix (cell × panel). Positive = case-aligned,
negative = control-aligned. `add_rfind_scores()` writes it into the Seurat metadata, one column per panel (like
`Seurat::AddModuleScore()`), so that the scores can be plotted and grouped like any other metadata.

From a Seurat object the raw counts are read from the `counts` layer of the `RNA` assay. If that layer is split
(e.g. after `merge()`), join it first: `mg <- SeuratObject::JoinLayers(mg, assay = "RNA")`.

### AUCell drop-in (single-list gene sets)

If you already use AUCell, RFind-sc accepts the same input shape via
`score_genesets()`:

```r
gene_sets <- list(
  microglia_marker = c("Tmem119", "P2ry12", "Cx3cr1"),
  inflammation     = c("Tnf", "Il1b", "Nfkb1")
)
out <- score_genesets(mg, gene_sets)   # cell x set scores
```

This is the AUCell-compatible entry point. To use RFind-sc's distinguishing
**direction-aware four-term score**, pass `(UP, DOWN)` panel pairs via
`score_panel_pairs()`:

```r
panels <- list(
  my_contrast = list(up = case_up_genes, down = case_down_genes)
)
out <- score_panel_pairs(mg, panels)   # signed scores
```

---

## Why RFind-sc

Existing per-cell scoring tools (AUCell, UCell, ssGSEA) treat gene sets as
**directionless markers** and return **non-negative** scores. That is
appropriate for pathway scoring against curated databases (MSigDB Hallmark,
GO), but limits analyses where gene sets come from case-vs-control contrasts —
the most common form of DEG output.

| | AUCell / UCell | RFind-sc |
|---|---|---|
| Uses UP and DOWN lists | ✗ (one list) | ✓ |
| Signed score (direction-aware) | ✗ (non-negative) | ✓ |
| Score | AUC / rank statistic | `-log10` of hypergeometric p (ranks only), cut-offs searched per cell |
| Baseline | per-cell ranks | expected counts from the dataset (P2 residual) |

**When to use which**: AUCell / UCell remain the right tools for directionless
pathway scoring (MSigDB Hallmark, GO). RFind-sc is the tool when your panels
have a case-vs-control structure (most DEG workflows).

---

## API

### Top-level

- `run_pipeline(input, biosets, config)` — Seurat or counts matrix + bioset list in, score matrix out.
- `score_genesets(input, gene_sets, config)` — **AUCell-compatible**: list of gene character vectors in, cell × set scores out.
  Each set is an unordered membership set and is used whole (gene order does not matter).
- `score_panel_pairs(input, panels, config)` — direction-aware: list of `(up, down)` gene character vectors in, signed cell × panel scores out.
  The lists are taken as ranked (most extreme first); add `ranked = FALSE` to a panel of unordered sets.
- `score_cells(counts, biosets, config, libsize, baseline_groups)` — lower-level (counts matrix only); same return
  shape. Pass `libsize` (total UMI over all genes) if `counts` has already been restricted to a gene universe, and
  `baseline_groups` (one label per cell, e.g. the cell type) for a per-group baseline (see below). `run_pipeline()`,
  `score_genesets()` and `score_panel_pairs()` accept both arguments too.
- `add_rfind_scores(object, scores, prefix = "RFind_")` — writes the scores of any of the functions above into the
  metadata of a Seurat object (or a data frame of cells), one column per panel, matched by cell name (0.3.7).
- `rfindsc_config(top_n, n_cell, nonzero_rate, clip, n_cores)` — hyperparameter defaults.

### Panel construction

- `build_panel_from_deg(deg, gene_col, log2fc_col, top_n, ...)` — DEG table → signed panel.
- `default_biosets()` — bundled 4 Hammond axes of the paper (DAM, Maturity, Interferon, G2M).
- `load_dedonno_biosets()` — bundled 4 De Donno axes of the paper (Glucocorticoid, Seizure, UPR, NF-kB).
- `msigdbr_to_panels(species, collection, subcollection)` — pull any MSigDB collection.
- `c2_cgp_auto_pair(species)` — pull C2 CGP `_UP/_DN` suffix sets and auto-pair them (~850 pairs).

### Core kernel

- `running_fisher(g1, g2, universe_size, clip, g2_ranked)` — Running Fisher score of two ranked lists: `-log10` of the
  smallest hypergeometric upper-tail p over every top-G × top-K cut-off pair. With `g2_ranked = FALSE`, `g2` is an
  unordered gene set: K is fixed at the set size and only the top-G cut-offs of `g1` are scanned, as for a biogroup in
  the original Running Fisher. (`bid_running_fisher()` is the former name, kept as an alias; before 0.3.0 it averaged
  two one-sided scans.)

---

## Panel tiers

Three ways to get panels into a scoring run:

1. **Bundled references of the paper**:
   ```r
   default_biosets()        # Hammond: DAM, Maturity, Interferon, G2M
   load_dedonno_biosets()   # De Donno: Glucocorticoid, Seizure, UPR, NF-kB
   ```
2. **Any MSigDB collection** on demand:
   ```r
   msigdbr_to_panels(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
   c2_cgp_auto_pair("Mus musculus")   # ~850 paired panels from C2 CGP
   ```
3. **Custom panel from your own DEG**:
   ```r
   my_panel <- build_panel_from_deg(
     deg, gene_col = "gene", log2fc_col = "log2FC",
     padj_col = "padj", padj_thresh = 0.05,
     label = "my_condition"
   )
   out <- run_pipeline(seurat_obj, biosets = list(my_condition = my_panel))
   ```
   The whole ranked table is kept; at scoring time the top `top_n` (500) per
   side are taken and then the genes of your dataset among them are used.

Combine tiers freely:
```r
biosets <- c(default_biosets(), load_dedonno_biosets(),
             list(my_condition = my_panel))
```

Each panel is scored in its pair universe: the genes of the dataset that are also in the panel's `eligible` genes (the genes
its source study reported; set by `build_panel_from_deg()` and carried by the bundled ranked references), or all genes of the
dataset for a panel without `eligible` (a gene set). Each cell's UP / DOWN lists are restricted to it before the top 500 are
taken (0.3.5); a ranked panel keeps its top 500 genes per side and then those in the pair universe (0.3.8), and a gene set
is used whole.

**Ranked lists and gene sets.** A panel is either ranked (DEG lists, most extreme first) or an unordered membership set
(`ranked = FALSE`: the bundled gene groups `Interferon`, `G2M`, `UPR`, `NF-kB`, every panel from
`msigdbr_to_panels()` / `c2_cgp_auto_pair()`, and every set given to `score_genesets()`). A set has no order, so it is
used whole: only the cut-offs of the cell's ranked lists are scanned (before 0.3.4 its alphabetical order was scanned
too). For a panel without a `ranked` field, types `"gene_group"`, `"msigdb_type_b"` and `"msigdb_cgp_paired"` count as
sets.

---

## Configuration

Defaults (see `?rfindsc_config`) follow the paper:

| parameter | default | meaning |
|---|---|---|
| `top_n` | 500 | cap per side (UP / DOWN) of each ranked panel, applied before the pair universe |
| `n_cell` | 500 | per-cell UP / DOWN signature length |
| `nonzero_rate` | 0.01 | gene universe filter (detected in ≥ 1% of cells) |
| `clip` | 300 | score magnitude cap (symmetric) |
| `n_cores` | cores − 1 (1 on Windows) | parallel mclapply workers |

Two environment variables tune memory and dispatch; scores do not depend on them. Unset, or set to a value that is not
a number, each takes its default:

| variable | default | meaning |
|---|---|---|
| `RFIND_DENSE_MAX_GB` | 8 | largest dense copy (GB; genes × cells × 8 bytes) of the counts made before the workers start; larger inputs are made dense one chunk at a time. Set it lower on a machine with little memory |
| `RFIND_CHUNK_SIZE` | 200 | largest number of cells per parallel task |

### Per-cell-type baseline

By default every cell is compared with one dataset-wide baseline (the mean CPM of all cells). With
`baseline_groups`, each cell's expected counts use the mean CPM of the cells of its own group instead, as the paper
does for the De Donno et al. hippocampus (22 cell types; both conditions pooled in each type's baseline), so that cell
identity does not dominate the signature:

```r
out <- run_pipeline(seu, biosets = load_dedonno_biosets(),
                    baseline_groups = seu$cell_type)   # one label per cell
```

The scoring universe is still chosen over all cells. Default `NULL` (one baseline; unchanged behaviour).

---

## Performance

About 3 s for 2,000 cells and the 4 bundled panels (`n_cell = 500`; 7 cores of an 8-core Apple M3 laptop; 8 s on one core; 2026-09-29).
Scales ≈ linearly in cell count.

---

## Companion tools

- **[rfindsc-py](https://github.com/tmurano/rfindsc-py)** — the Python implementation of this package (same
  version numbers, scores identical to RFindsc to < 1e-8; tested on shared fixtures).
- **[RFind-sc](https://github.com/tmurano/RFind-sc)** — the analysis code of the paper, including the scripts that build
  the bundled references (`analysis/hammond/06_build_references.R`, `analysis/dedonno/03_references.R`).

---

## Citation

Murano et al. 2026 (in preparation). Please cite the preprint once available.
For now, reference this repository and the version used (the paper used 0.3.8, git tag `v0.3.8`).

---

## License

MIT © 2026 Tomoyuki Murano (Miyakawa Lab)
