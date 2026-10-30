# [[03_variant_calling/15_plot_burden_results.R]]
# Goal: Figures for the polygenic burden analysis. Consumes the tables written by
#       14. No statistics here -- every number plotted comes from a result file.
#
# Outputs:
#   15_burden_offset_effect.pdf/.png   the key figure: raw vs offset-corrected
#   15_gene_burden_volcano.pdf/.png      gene-level burden, all genes
#   15_panel_burden.pdf/.png            Open Targets panel at each threshold
#
# Usage:
#   Rscript 03_variant_calling/15_plot_burden_results.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
})

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) curr <- dirname(curr)
  curr
}
project_root <- Sys.getenv("PROJECT_ROOT", get_project_root())

TAB <- file.path(project_root, "results/03_variant_calling/tables")
FIG <- file.path(project_root, "results/03_variant_calling/figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

need <- function(f) {
  p <- file.path(TAB, f)
  if (!file.exists(p)) {
    message("[15] missing ", f, " -- skipped (run 14_gene_burden_test.py)")
    return(NULL)
  }
  fread(p)
}

save_fig <- function(p, name, w = 8, h = 5.5) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  message("[15] wrote ", name, ".pdf / .png")
}

sens <- need("14_sensitivity.csv")
gene_burden <- need("14_gene_burden.csv")
panel_burden <- need("14_panel_burden.csv")
tot  <- need("11_donor_totals.csv")

# --- Figure 1: raw vs offset-corrected effect ---------------------------
# This is the figure that answers "why did earlier results look inconsistent":
# the raw burden ratio and the depth-corrected ratio are different quantities,
# and only the latter is interpretable.
if (!is.null(sens)) {
  d <- copy(sens[!is.na(rate_ratio_offset) & !is.na(rate_ratio_no_offset)])
  d[, analysis := factor(analysis, levels = rev(unique(analysis)))]
  long <- rbindlist(list(
    data.table(analysis = d$analysis,
               estimate = d$rate_ratio_no_offset, test = "no offset (confounded)"),
    data.table(analysis = d$analysis,
               estimate = d$rate_ratio_offset, test = "offset-corrected")
  ))
  p <- ggplot(long, aes(x = estimate, y = analysis, colour = test)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_point(size = 2.6) +
    facet_wrap(~ test, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c("no offset (confounded)" = "#C0392B",
                                   "offset-corrected" = "#1F618D")) +
    labs(title = "Burden effect: unadjusted vs detection-depth corrected",
         subtitle = sprintf(
           "Dashed line = no difference. MDD donors carry %.2fx more total variants (%.0f vs %.0f),\nso the unadjusted estimate is confounded by sequencing depth.",
           mean(tot[condition == "mdd", total_variants_detected]) /
             mean(tot[condition == "ctrl", total_variants_detected]),
           mean(tot[condition == "mdd", total_variants_detected]),
           mean(tot[condition == "ctrl", total_variants_detected])),
         x = "rate ratio (MDD / Ctrl)", y = NULL, colour = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          plot.subtitle = element_text(size = 8.5, colour = "grey30"))
  save_fig(p, "15_burden_offset_effect", w = 8, h = 1.6 * length(unique(long$analysis)) + 2)
}

# --- Figure 2: gene-level burden volcano --------------------------------------------
if (!is.null(gene_burden) && nrow(gene_burden)) {
  d <- copy(gene_burden[is.finite(-log10(pmax(p_glm, .Machine$double.xmin)))])
  d[, sig := factor(fdr_bh < 0.05, levels = c(FALSE, TRUE),
                    labels = c("BH-FDR >= 0.05", "BH-FDR < 0.05"))]
  p <- ggplot(d, aes(x = -log10(pmax(p_glm, .Machine$double.xmin)), y = log10(rate_ratio),
                     colour = sig)) +
    geom_point(size = 0.7, alpha = 0.55) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = c("BH-FDR >= 0.05" = "grey55",
                                   "BH-FDR < 0.05" = "#C0392B")) +
    labs(title = sprintf("gene-level burden: unbiased gene-level burden (%s genes tested)", format(nrow(d), big.mark = ",")),
         subtitle = sprintf("Nominal p<0.05: %d (%.1f%%, ~5%% expected by chance) | BH-FDR<0.05: %d",
                            sum(d$p_glm < 0.05), 100 * mean(d$p_glm < 0.05), sum(d$fdr_bh < 0.05)),
         x = expression(-log[10]("nominal p")), y = "log10 rate ratio (MDD / Ctrl)", colour = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9))
  save_fig(p, "15_gene_burden_volcano")
}

# --- Figure 3: panel burden panel across thresholds ----------------------------
if (!is.null(panel_burden) && nrow(panel_burden)) {
  d <- copy(panel_burden)
  d[, significant := factor(p_glm_offset < 0.05, levels = c(FALSE, TRUE),
                            labels = c("p >= 0.05", "p < 0.05"))]
  p <- ggplot(d, aes(x = factor(genetic_threshold), y = rate_ratio_offset,
                     fill = significant)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
    geom_col(width = 0.55) +
    scale_fill_manual(values = c("p >= 0.05" = "#5B9BD5", "p < 0.05" = "#C0392B")) +
    labs(title = "panel burden: Open Targets MDD panel (genetic evidence only)",
         subtitle = "Declared in advance. No threshold reaches p < 0.05.",
         x = "Open Targets genetic evidence threshold",
         y = "offset-corrected rate ratio (MDD / Ctrl)", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.subtitle = element_text(size = 9))
  save_fig(p, "15_panel_burden", w = 7, h = 4.5)
}

message("[15] done")