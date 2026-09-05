# [[03_variant_calling/13_asm_compare_viz.R]]
# Goal: Produce ASM visualisation figures from the compute artifact written by 12_asm_compare.R (asm_results.rds).

suppressPackageStartupMessages({
  suppressMessages({
    ok_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  })
})
if (!ok_ggplot) stop("ggplot2 is required.")

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

rds_path <- file.path(project_root, "_data/RData/03_variant_calling/asm_results.rds")
if (!file.exists(rds_path)) {
  stop(paste0("Artifact not found: ", rds_path, ". Run 12_asm_compare.R first."))
}

cfg <- readRDS(rds_path)
res <- cfg$results
sweep <- cfg$sensitivity
P <- cfg$params

fig_dir <- file.path(project_root, "results/03_variant_calling/figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

tag       <- ifelse(isTRUE(P$input_available), "Input-null offset", "0.5-null offset")
xlabel    <- "coef (log-OR MDD/control)"
ylabel    <- "-log10(BH-adjusted p)"

# --- 1. Volcano Plot ---
res$coef_log2 <- res$coef / log(2)             # natural-log OR -> log2
res$fdr_neglog <- -log10(res$fdr)              # NA stays NA
sig <- !is.na(res$fdr) & res$fdr < P$fdr_cut
res$direction <- ifelse(res$coef > 0, "up", "down")
res$group <- ifelse(sig, paste0("sig-", res$direction), "n.s.")

vol <- ggplot2::ggplot(res, ggplot2::aes(x = .data$coef_log2, y = .data$fdr_neglog,
                                         color = .data$group))
vol <- vol + ggplot2::geom_point(alpha = 0.7, size = 1.6, na.rm = TRUE)
vol <- vol + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40")
sig_txt <- res[sig, , drop = FALSE]
has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

if (has_ggrepel) {
  vol <- vol + ggrepel::geom_text_repel(
    data = sig_txt,
    ggplot2::aes(label = .data$gene),
    size = 3.6,
    fontface = "bold.italic",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.4,
    max.overlaps = Inf,
    show.legend = FALSE,
    na.rm = TRUE
  )
} else {
  # Fallback: Manual offset nudging if ggrepel is absent
  sig_txt$nudge_y <- 0.05
  sig_txt$nudge_x <- 0.0
  sig_txt$nudge_y[sig_txt$gene == "BIN3"] <- 0.12
  sig_txt$nudge_y[sig_txt$gene == "MGAT5"] <- -0.08
  sig_txt$nudge_x[sig_txt$gene == "BIN3"] <- -0.15
  sig_txt$nudge_x[sig_txt$gene == "MGAT5"] <- 0.15
  
  vol <- vol + ggplot2::geom_text(
    data = sig_txt,
    ggplot2::aes(label = .data$gene, x = .data$coef_log2 + .data$nudge_x, y = .data$fdr_neglog + .data$nudge_y),
    size = 3.4, fontface = "bold.italic", show.legend = FALSE, na.rm = TRUE)
}
vol <- vol + ggplot2::scale_color_manual(
  "MDD vs Ctrl",
  values = c("sig-up" = "#c0392b", "sig-down" = "#2980b9", "n.s." = "grey70"),
  breaks = c("sig-up", "sig-down", "n.s."),
  labels = c("sig-up" = "up (MDD>Ctrl)", "sig-down" = "down (Ctrl>MDD)", "n.s." = "not significant"))
vol <- vol + ggplot2::theme_bw(base_size = 12)
vol <- vol + ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal",
                            legend.key.size = grid::unit(0.9, "lines"))
vol <- vol + ggplot2::labs(
  x = xlabel, y = ylabel,
  title = sprintf("ASM loci: MDD vs Control (%s)", tag),
  caption = sprintf("BH-FDR cutoff = %.2f, SKEW_CUT = %.1f, min subjects = %d",
                    P$fdr_cut, P$skew_cut, P$min_subjects))

out_pdf <- file.path(fig_dir, "13_asm_volcano.pdf")
out_png <- file.path(fig_dir, "13_asm_volcano.png")
ggplot2::ggsave(out_pdf, vol, width = 180, height = 140, units = "mm", dpi = 300)
ggplot2::ggsave(out_png, vol, width = 180, height = 140, units = "mm", dpi = 300)
message("[INFO] Wrote volcano plot to: ", out_pdf, " and ", out_png)

# --- 2. Sensitivity Sweep Plot (Barplot) ---
if (nrow(sweep) > 0) {
  sweep$skew_lbl <- sprintf("skew >= %.1f", sweep$skew_cut)
  sweep$skew_lbl <- factor(sweep$skew_lbl,
                           levels = sprintf("skew >= %.1f", sort(unique(sweep$skew_cut))))
  sp <- ggplot2::ggplot(sweep, ggplot2::aes(x = factor(.data$min_subjects), y = .data$n_loci,
                                            fill = .data$skew_lbl))
  sp <- sp + ggplot2::geom_col(position = "dodge", alpha = 0.9)
  sp <- sp + ggplot2::geom_text(
    ggplot2::aes(label = .data$n_loci, group = .data$skew_lbl),
    position = ggplot2::position_dodge(width = 0.9), size = 3,
    vjust = -0.5)
  sp <- sp + ggplot2::theme_bw(base_size = 12)
  sp <- sp + ggplot2::scale_fill_brewer(palette = "Set2")
  sp <- sp + ggplot2::expand_limits(y = max(sweep$n_loci) * 1.1)
  sp <- sp + ggplot2::labs(x = "min subjects", y = "tested loci", fill = NULL,
                           title = "ASM Sensitivity to Filtering Parameters")
  
  out_sp_pdf <- file.path(fig_dir, "13_asm_sensitivity.pdf")
  out_sp_png <- file.path(fig_dir, "13_asm_sensitivity.png")
  ggplot2::ggsave(out_sp_pdf, sp, width = 180, height = 140, units = "mm", dpi = 300)
  ggplot2::ggsave(out_sp_png, sp, width = 180, height = 140, units = "mm", dpi = 300)
  message("[INFO] Wrote sensitivity plot to: ", out_sp_pdf, " and ", out_sp_png)
} else {
  message("[WARN] sensitivity table empty; skipping sweep plot.")
}

cat("[INFO] Done. ASM visualization figures saved in ", fig_dir, "\n", sep = "")
