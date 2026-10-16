#' Convert an MSigDB collection to RFind-sc bioset form (Type B)
#'
#' Pulls any MSigDB collection via `msigdbr` and wraps each gene set as a Type B
#' panel (UP only; DOWN empty). Useful for pathway scoring against the full
#' MSigDB library. MSigDB gene sets are membership sets with no order, so each
#' panel carries `ranked = FALSE` and is used whole at scoring (see
#' [running_fisher()]; 0.3.4).
#'
#' @param species Character: `msigdbr` species name (e.g. `"Mus musculus"`,
#'   `"Homo sapiens"`). Default `"Mus musculus"`.
#' @param collection Character: MSigDB collection code. `"H"` for Hallmark,
#'   `"C2"` for curated (includes CGP / CP), `"C5"` for GO, etc.
#' @param subcollection Optional character: sub-collection code (e.g. `"CGP"`
#'   within `"C2"`).
#' @param strip_prefix Logical: strip the collection prefix (e.g. `"HALLMARK_"`)
#'   from panel names. Default TRUE.
#'
#' @return A named list of Type B biosets (with `ranked = FALSE`). Requires
#'   `msigdbr` (>= 10.0.0, for the `collection` / `subcollection` arguments)
#'   installed.
#' @export
msigdbr_to_panels <- function(species = "Mus musculus",
                              collection = "H",
                              subcollection = NULL,
                              strip_prefix = TRUE) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required for msigdbr_to_panels().")
  }
  args <- list(species = species, collection = collection)
  if (!is.null(subcollection)) args$subcollection <- subcollection
  df <- do.call(msigdbr::msigdbr, args)
  gs_names <- sort(unique(df$gs_name))
  out <- lapply(gs_names, function(gs) {
    genes <- unique(df$gene_symbol[df$gs_name == gs])
    list(up = data.frame(gene = genes,
                         diff = rep(NA_real_, length(genes)),
                         rank = seq_along(genes),
                         stringsAsFactors = FALSE),
         down = empty_panel_side(),
         type = "msigdb_type_b",
         ranked = FALSE,
         label = gs,
         source = sprintf("MSigDB %s%s (%s)",
                          collection,
                          if (is.null(subcollection)) "" else paste0(":", subcollection),
                          species))
  })
  panel_names <- if (isTRUE(strip_prefix)) {
    sub(sprintf("^%s_", toupper(collection)),
        "",
        sub("^HALLMARK_", "", gs_names))
  } else gs_names
  names(out) <- panel_names
  out$meta <- list(source = sprintf("msigdbr %s / %s%s", species, collection,
                                    if (is.null(subcollection)) ""
                                    else paste0("/", subcollection)),
                   n_panels = length(gs_names),
                   fetched_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  out
}

#' Auto-pair MSigDB C2 CGP _UP / _DN gene sets into bidirectional panels
#'
#' The MSigDB C2 Chemical & Genetic Perturbations (CGP) collection contains
#' many gene sets with paired `_UP` / `_DN` suffixes (the two sides of a single
#' perturbation DEG). This helper pulls C2 CGP and auto-pairs them into Type A
#' (bidirectional) panels: roughly 850 pairs (about 1,700 individual gene sets)
#' as of MSigDB 2024. The MSigDB sets carry no gene order, so each panel has
#' `ranked = FALSE` and both sides are used whole at scoring (0.3.4).
#'
#' @param species Character. Default `"Mus musculus"`.
#' @param up_suffix,dn_suffix Character suffixes to strip to obtain the panel
#'   stem. Defaults `"_UP"` and `"_DN"`.
#'
#' @return A named list of Type A biosets (each with populated `$up` and
#'   `$down`, and `ranked = FALSE`).
#' @export
c2_cgp_auto_pair <- function(species = "Mus musculus",
                             up_suffix = "_UP",
                             dn_suffix = "_DN") {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required for c2_cgp_auto_pair().")
  }
  df <- msigdbr::msigdbr(species = species,
                         collection = "C2",
                         subcollection = "CGP")
  gs_names <- unique(df$gs_name)
  up_pat <- paste0(up_suffix, "$")
  dn_pat <- paste0(dn_suffix, "$")
  up_set <- grep(up_pat, gs_names, value = TRUE)
  dn_set <- grep(dn_pat, gs_names, value = TRUE)
  up_stem <- sub(up_pat, "", up_set)
  dn_stem <- sub(dn_pat, "", dn_set)
  shared <- intersect(up_stem, dn_stem)

  biosets <- lapply(shared, function(stem) {
    up_genes <- unique(df$gene_symbol[df$gs_name == paste0(stem, up_suffix)])
    dn_genes <- unique(df$gene_symbol[df$gs_name == paste0(stem, dn_suffix)])
    list(up = data.frame(gene = up_genes,
                         diff = rep(NA_real_, length(up_genes)),
                         rank = seq_along(up_genes),
                         stringsAsFactors = FALSE),
         down = data.frame(gene = dn_genes,
                           diff = rep(NA_real_, length(dn_genes)),
                           rank = seq_along(dn_genes),
                           stringsAsFactors = FALSE),
         type = "msigdb_cgp_paired",
         ranked = FALSE,
         label = stem,
         source = sprintf("MSigDB C2:CGP auto-paired (%s / %s)",
                          up_suffix, dn_suffix))
  })
  names(biosets) <- shared
  biosets$meta <- list(
    source = sprintf("msigdbr %s / C2:CGP auto-paired", species),
    n_panels = length(shared),
    n_unpaired = length(gs_names) - 2 * length(shared),
    fetched_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  )
  biosets
}

#' Empty panel side (helper)
#' @keywords internal
empty_panel_side <- function() {
  data.frame(gene = character(0), diff = numeric(0), rank = integer(0),
             stringsAsFactors = FALSE)
}
