# [[03_variant_calling/21_plot_burden_control_vs_mdd.R]]
# Goal: Shared plotting helper for the per-collection burden figures
#       (22_plot_burden_reactome.R, 23_*_kegg.R, 24_*_go_bp.R, 25_*_go_cc.R,
#        26_*_go_mf.R). NOT run directly -- sourced by those scripts.
#
# Every figure has the same two-panel layout:
#   LEFT   the ten most impacted sets OF THIS COLLECTION. Scoped to the
#          collection so the two panels describe one database and its BH-FDR
#          denominator applies to both, and so no half of the figure is a copy
#          shared with the other four.
#   RIGHT  the synapse / neurotransmission sets of that same collection, with
#          the collection's own test count stated.
# Panel A ordering by |effect| is the honest way to present a null: it shows
# what the extremes of noise look like.
#
# Styling follows utils/plot_gtex_gwas_integration.R and plot_gtex_eqtl_summary.R:
#   theme_minimal(base_family="sans"), slate palette, bold titles, bg="white".
#
# Inputs (written by 15, 16b, 17, 20):
#   20_burden_per_donor.csv
#   20_burden_summary.csv

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

PAL_CTRL <- "#1e3a8a"   # navy
PAL_MDD  <- "#dc2626"   # red
INK      <- "#0f172a"
INK_SOFT <- "#64748b"
PANEL_BG <- "#fafbfc"
GRID     <- "#e2e8f0"

# Minimum genes a set must overlap before it is eligible for the LEFT panel.
# Without this floor the top-10 is filled by sets overlapping 3 genes, where the
# burden ratio is a quotient of two or three numbers and carries no information
# -- it is denominator noise dressed up as an effect. 10 overlapping genes is
# still a weak set, but it is the smallest denominator at which a group
# comparison means anything, and it is the same floor for every collection.
MIN_GENES_PANEL <- 10

# Largest a displayed panel may be. GO_BP has 342 sets matching the neuro
# pattern; drawing all of them produces an unreadable strip.
MAX_PANEL_ROWS  <- 14

# Below this many matching sets a collection panel is a near-empty strip rather
# than a comparison -- KEGG matches exactly 1 set on the neuro vocabulary, which
# reads as a broken panel. Under this floor the panel falls back to the
# collection's largest effects and the subtitle says so.
MIN_NEURO_PANEL <- 5

pretty_name <- function(x) {
  x <- sub("^(HALLMARK|REACTOME|KEGG|WP|BIOCARTA|GOBP|GOCC|GOMF)_", "", x)
  x <- gsub("_", " ", x)
  # strsplit on "_" leaves fragments like "O_GLYCAN" whose second word starts
  # lowercase; toTitleCase would render them as "o Glycan". Title-case each word
  # and restore the acronyms that gene-family names are conventionally written in.
  vapply(strsplit(x, " "), function(w) {
    w <- tolower(w)
    w <- paste0(toupper(substr(w, 1, 1)), substring(w, 2))
    keep <- c("Gaba", "Nmda", "Ampa", "Gq", "Gabaa", "Nmda", "Cns", "Dna", "Rna")
    w[tolower(w) %in% tolower(keep)] <- toupper(w[tolower(w) %in% tolower(keep)])
    paste(w, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}

wrap_label <- function(x, width = 32) {
  # Truncate before wrapping, not after. strwrap() keeps every character, so a
  # very long GO term still needs four lines and collides with the next row's
  # label in a fixed-height strip. 60 chars is enough to identify any GO term
  # (they are named, not numbered, in the table).
  x <- ifelse(nchar(x) > 60, paste0(substr(x, 1, 57), "..."), x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

# Synapse / neurotransmission vocabulary per collection. Patterns, not explicit
# lists, so a collection that has no such set is handled rather than crashing.
NEURO_PAT <- c(
  HALLMARK    = "SYNAP|NEURON|NEUROTRANSMIT|GLUTAMAT|GABA",
  REACTOME    = "SYNAP|NEURON|NEUROTRANSMIT|GLUTAMAT|GABA|NMDA|KAINATE|POSTSYNAPTIC|DOPAMIN",
  KEGG        = "GLUTAMAT|GABA|SYNAP|NEURON",
  WIKIPATHWAYS = "SYNAP|NEURON|GLUTAMAT|GABA|NEUROTRANSMIT",
  BIOCARTA    = "SYNAP|NEURON|GABA",
  GO_BP       = "SYNAP|NEURON|NEUROTRANSMIT|GLUTAMAT|GABA|AXON|DENDRIT",
  GO_CC       = "SYNAP|NEURON|DENDRIT|AXON|NEURO",
  GO_MF       = paste("SYNAP|NEUROTRANSMIT|GLUTAMAT|GABA|NMDA|KAINATE|CHOLINERGIC",
                      "|ADRENERGIC|DOPAMINERGIC|SEROTONINERGIC", collapse = "")
)

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) curr <- dirname(curr)
  curr
}

load_burden <- function() {
  root <- get_project_root()
  ctx <- list(TAB = file.path(root, "results/03_variant_calling/tables"),
              FIG = file.path(root, "results/03_variant_calling/figures"))
  pd <- fread(file.path(ctx$TAB, "20_burden_per_donor.csv"))
  sm <- fread(file.path(ctx$TAB, "20_burden_summary.csv"))
  pd[, grp := factor(ifelse(condition == "mdd", "MDD (n=20)", "Control (n=19)"),
                     levels = c("Control (n=19)", "MDD (n=20)"))]
  pd[, label := wrap_label(pretty_name(set_name), 32)]
  sm[, label := wrap_label(pretty_name(set_name), 32)]
  list(pd = pd, sm = sm, ctx = ctx)
}

# Shared strip geometry so the LEFT and RIGHT panels line up pixel-for-pixel.
strip_plot <- function(d, title, subtitle, x_title) {
  ggplot(d, aes(x = burden, y = label)) +
    geom_vline(xintercept = 0, color = GRID, linewidth = 0.6) +
    geom_point(aes(color = grp), size = 1.7, alpha = 0.75,
               position = position_jitter(height = 0.13, width = 0)) +
    stat_summary(aes(group = grp), fun = median, geom = "point", size = 2.8,
                 shape = 18, color = INK, position = position_dodge(width = 0.55)) +
    scale_color_manual(values = c("Control (n=19)" = PAL_CTRL, "MDD (n=20)" = PAL_MDD)) +
    labs(title = title, subtitle = strwrap(subtitle, width = 52),
         x = x_title, y = NULL) +
    theme_minimal(base_family = "sans") +
    theme(
      plot.title = element_text(size = 15, face = "bold", color = INK, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 11.5, face = "italic", color = INK_SOFT, margin = margin(b = 12)),
      axis.text.y = element_text(size = 10.5, color = INK, lineheight = 0.92),
      axis.text.x = element_text(size = 11.5, color = INK),
      axis.title.x = element_text(size = 12.5, face = "bold", color = INK, margin = margin(t = 8)),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = GRID, linewidth = 0.5, linetype = "dotted"),
      panel.background = element_rect(fill = PANEL_BG, color = "#cbd5e1", linewidth = 0.7),
      legend.position = "none",
      plot.margin = margin(t = 10, r = 14, b = 10, l = 8)
    )
}

strip_of <- function(pd, sel) {
  # No merge: pd already carries the label derived from set_name, and merging
  # sel's copy back in would produce label.x / label.y for no gain.
  d <- pd[set_name %chin% sel$set_name]
  d[, label := factor(label, levels = rev(sel$label))]
  d
}

top_panel <- function(pd, sm, collection) {
  coll <- data.table::copy(collection)   # stop data.table scoping the argument
  sel <- sm[collection == coll & genes_in_cohort >= MIN_GENES_PANEL]
  sel <- sel[order(-abs(log2rr))][1:10]
  strip_plot(strip_of(pd, sel),
             sprintf("A   %s: ten most impacted gene sets", coll),
             sprintf("largest effects in this collection (>= %d genes overlapping)",
                     MIN_GENES_PANEL),
             "variants carried")
}

collection_panel <- function(pd, sm, collection, title, subtitle) {
  coll <- data.table::copy(collection)   # stop data.table scoping the argument
  in_coll <- sm[collection == coll]
  n_tested <- nrow(in_coll)

  pat <- NEURO_PAT[[coll]]
  sel <- in_coll[grepl(pat, set_name)]
  if (nrow(sel) < MIN_NEURO_PANEL) {
    # HALLMARK has no synapse or neurotransmission set at all, and KEGG matches
    # only one. Fall back to the largest effects the collection does have, and
    # say so in the subtitle rather than drawing a one-row strip.
    # Rows 11+ of the ranking: rows 1-10 are already Panel A, and repeating them
    # here would make the two halves of the same figure identical.
    sel <- in_coll[genes_in_cohort >= MIN_GENES_PANEL][order(-abs(log2rr))]
    sel <- sel[-(1:10)][1:MAX_PANEL_ROWS]
    # Rename the panel: the fallback is not synapse biology, and a title that
    # claims otherwise would misrepresent what is drawn.
    title <- sprintf("B   %s: ranks 11-%d by |effect|", coll, 10 + nrow(sel))
    subtitle <- sprintf("%s -- only %d set in this collection matches the synapse vocabulary",
                        subtitle, nrow(in_coll[grepl(pat, set_name)]))
  } else if (nrow(sel) > MAX_PANEL_ROWS) {
    # Keep the biggest absolute effects; all displayed sets are in the table.
    sel <- sel[genes_in_cohort >= 3][order(-abs(log2rr))][1:MAX_PANEL_ROWS]
  } else {
    sel <- sel[order(log2rr)]
  }
  if (!nrow(sel)) stop("no sets to plot for collection: ", coll)

  d <- strip_of(pd, sel)
  strip_plot(d, title, subtitle,
             sprintf("variants carried   (%s gene sets tested)", format(n_tested, big.mark = ",")))
}

compose_figure <- function(pA, pB) {
  p <- pA + pB + plot_layout(guides = "collect") &
    theme(legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 12.5, face = "bold", color = "#334155"))
  p + plot_annotation(
    caption = paste0(
      "Each point is one donor (39 total); diamonds are group medians.   ",
      "Control and MDD overlap in every set shown.   ",
      "BH-FDR < 0.05: 0 gene sets significant in any collection.   ",
      "EXPLORATORY analysis -- supports no claim on its own."),
    theme = theme(
      plot.caption = element_text(size = 10.5, face = "italic", color = INK_SOFT,
                                  hjust = 0, margin = margin(t = 12))))
}

save_figure <- function(p, ctx, name, n_rows) {
  dir.create(ctx$FIG, recursive = TRUE, showWarnings = FALSE)
  h <- max(5.0, 0.40 * n_rows + 3.4)
  ggsave(file.path(ctx$FIG, paste0(name, ".pdf")), p,
         width = 15.5, height = h, dpi = 300, bg = "white")
  ggsave(file.path(ctx$FIG, paste0(name, ".png")), p,
         width = 15.5, height = h, dpi = 300, bg = "white")
  message("[", name, "] wrote ", name, ".pdf / .png")
}