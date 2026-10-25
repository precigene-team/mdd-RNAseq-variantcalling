# [[02_dge_analysis/13_gsea_figures.R]]
# Goal: Publication figures for the preranked Hallmark GSEA run in 12_gsea_fgsea.R.
#       Aesthetic follows 08_enrichment_plotting.R (red/blue NES bars, theme_bw).
#
# Steps:
#   1. Load 12_<contrast>_gsea_hallmark.csv for the two input contrasts.
#   2. Fig A: NES barplot per contrast (top 12 by p-value, red = positive NES,
#      blue = negative NES), one file per contrast.
#   3. Fig B: merged dotplot - both cohorts on the same set rows, with a
#      dumbbell segment joining each set's male and female NES, so agreement
#      and overlap between cohorts read at a glance.
#   4. Fig C: fgsea enrichment plots for the top 3 sets in the male contrast,
#      with leading-edge shading, as a 3-panel patchwork.
#   5. Fig D: male-vs-female NES scatter over all 50 sets - shows that the
#      female contrast sits on the null while the male contrast does not.
#   6. Save PNGs to results/02_dge_analysis/figures/ prefixed 13_.

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(stringr)
  library(patchwork)
  library(fgsea)
  library(scales)
  library(tidyr)
})

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

tables_dir <- file.path(project_root, "results/02_dge_analysis/tables")
gsea_dir   <- file.path(project_root, "results/02_dge_analysis/gsea")
figures_dir <- file.path(project_root, "results/02_dge_analysis/figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

clean_pathway_name <- function(names_vec) {
  x <- gsub("^HALLMARK_", "", names_vec)
  x <- gsub("_", " ", x)
  x <- stringr::str_to_title(x)
  x <- gsub("Nfkb", "NFkB", x)
  x <- gsub("Il6", "IL6", x)
  x <- gsub("Jak", "JAK", x)
  x <- gsub("Stat3", "STAT3", x)
  x <- gsub("Stat5", "STAT5", x)
  x <- gsub("Dna", "DNA", x)
  x
}

# Same palette as 08_enrichment_plotting.R
pal <- c("Positive NES" = "#e31a1c", "Negative NES" = "#1f78b4")

base_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9, face = "bold", lineheight = 0.85),
        axis.title.y = element_blank(),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.position = "bottom")

load_gsea <- function(contrast) {
  f <- file.path(tables_dir, paste0("12_", contrast, "_gsea_hallmark.csv"))
  if (!file.exists(f)) stop("missing table: ", f)
  read_csv(f, show_col_types = FALSE) %>%
    mutate(pathway_clean = clean_pathway_name(pathway),
           neg_log_p = -log10(pmax(pvalue, 1e-10)),
           Direction = ifelse(NES > 0, "Positive NES", "Negative NES"))
}

# -----------------------------------------------------------------------------
# Fig A - NES barplot per contrast
# -----------------------------------------------------------------------------
nes_barplot <- function(df, contrast, top_n = 12) {
  plot_df <- df %>%
    mutate(pathway_clean = str_wrap(pathway_clean, 34)) %>%
    arrange(pvalue) %>%
    head(top_n)

  p <- ggplot(plot_df, aes(x = NES, y = reorder(pathway_clean, NES),
                           fill = Direction)) +
    geom_col(width = 0.75, alpha = 0.85) +
    scale_fill_manual(values = pal, name = NULL) +
    geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
    geom_text(aes(label = sprintf("%.2f", NES)),
              hjust = ifelse(plot_df$NES > 0, -0.15, 1.15),
              size = 3.0, color = "grey15") +
    base_theme +
    labs(title = paste0("Hallmark GSEA - ", gsub("^deg_", "", gsub("_", " ", contrast))),
         subtitle = "Preranked, Wald statistic ranking, MSigDB Hallmark (50 sets)",
         x = "Normalized Enrichment Score (NES)")

  height <- max(4.5, nrow(plot_df) * 0.42)
  ggsave(file.path(figures_dir, paste0("13_", contrast, "_gsea_hallmark_nes.png")),
         plot = p, width = 8.5, height = height, dpi = 300, bg = "white",
         limitsize = FALSE)
  message("  wrote 13_", contrast, "_gsea_hallmark_nes.png")
}

# -----------------------------------------------------------------------------
# Fig B - merged dotplot: both cohorts on the same set rows
# -----------------------------------------------------------------------------
male_df   <- load_gsea("deg_male_input")
female_df <- load_gsea("deg_female_input")

nes_barplot(male_df, "deg_male_input")
nes_barplot(female_df, "deg_female_input")

dots <- bind_rows(
  mutate(male_df,   contrast = "Male Input",   neg_log_padj = -log10(padj)),
  mutate(female_df, contrast = "Female Input", neg_log_padj = -log10(padj))
)

# rows = all 50 sets, ordered by the male NES so the male gradient reads top to
# bottom; alpha marks which sets cleared padj < 0.05
lev <- male_df %>% arrange(NES) %>% pull(pathway_clean)

dots <- dots %>%
  mutate(pathway_clean = str_wrap(pathway_clean, 26),
         pathway_clean = factor(pathway_clean, levels = str_wrap(lev, 26)),
         sig = padj < 0.05)

# one segment per set, joining its male and female NES
segs <- dots %>%
  select(pathway_clean, contrast, NES) %>%
  pivot_wider(names_from = contrast, values_from = NES) %>%
  filter(!is.na(`Male Input`), !is.na(`Female Input`)) %>%
  transmute(pathway_clean, x = `Male Input`, xend = `Female Input`,
            y = pathway_clean, yend = pathway_clean)

p_dot <- ggplot(dots, aes(x = NES, y = pathway_clean, colour = contrast,
                          size = neg_log_padj, alpha = sig)) +
  geom_segment(data = segs, aes(x = x, xend = xend, y = y, yend = yend),
               inherit.aes = FALSE, colour = "grey70", linewidth = 0.9) +
  geom_point(shape = 21, position = position_dodge(width = 0.45)) +
  scale_colour_manual(values = c("Male Input" = "#2673A6",
                                 "Female Input" = "#D9592E"), name = NULL) +
  scale_size(range = c(2.5, 6), name = "-log10(padj)") +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey30") +
  base_theme +
  theme(axis.text.y = element_text(size = 8, face = "bold", lineheight = 0.8),
        axis.title.y = element_blank(),
        legend.position = "right") +
  labs(title = "Hallmark GSEA - Male and Female Input Cohorts",
       subtitle = "All 50 sets, joined by set; solid point = padj < 0.05",
       x = "Normalized Enrichment Score (NES)")

ggsave(file.path(figures_dir, "13_gsea_hallmark_dotplot.png"),
       plot = p_dot, width = 10.5, height = 7.2, dpi = 300, bg = "white")
message("  wrote 13_gsea_hallmark_dotplot.png")

# -----------------------------------------------------------------------------
# Fig C - enrichment plots for the top 3 male sets
# -----------------------------------------------------------------------------
stats_male <- {
  lines <- readLines(file.path(gsea_dir, "deg_male_input.rnk"))
  parts <- strsplit(lines[!startsWith(lines, "#")], "\t")
  s <- as.numeric(vapply(parts, `[`, character(1), 2))
  names(s) <- vapply(parts, `[`, character(1), 1)
  s[order(-s)]
}
gene_sets <- gmtPathways(file.path(gsea_dir, "hallmark.gmt"))

top3 <- male_df %>% arrange(pvalue) %>% head(3)

plots <- lapply(top3$pathway, function(pw) {
  plotEnrichment(gene_sets[[pw]], stats_male) +
    scale_y_continuous(name = "Running Enrichment Score (NES)") +
    scale_x_continuous(name = "Ranked genes (Wald statistic)") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 9, hjust = 0.5),
          legend.position = "none") +
    labs(title = clean_pathway_name(pw),
         subtitle = sprintf("NES = %.2f, padj = %.2e",
                            male_df$NES[match(pw, male_df$pathway)],
                            male_df$padj[match(pw, male_df$pathway)]))
})

p_enrich <- plots[[1]] + plots[[2]] + plots[[3]] +
  plot_layout(nrow = 1) +
  plot_annotation(title = "Leading-Edge Enrichment Plots - Male Input Cohort",
                  theme = theme(plot.title = element_text(size = 13, face = "bold",
                                                          hjust = 0.5)))

ggsave(file.path(figures_dir, "13_deg_male_input_gsea_hallmark_enrichment.png"),
       plot = p_enrich, width = 11.5, height = 4.0, dpi = 300, bg = "white")
message("  wrote 13_deg_male_input_gsea_hallmark_enrichment.png")

# -----------------------------------------------------------------------------
# Fig D - male vs female NES scatter over all 50 sets
# -----------------------------------------------------------------------------
scatter_df <- male_df %>%
  select(pathway, pathway_clean, NES_male = NES, padj_male = padj) %>%
  inner_join(female_df %>% select(pathway, NES_female = NES, padj_female = padj),
            by = "pathway") %>%
  mutate(pathway_clean = str_wrap(pathway_clean, 22),
           sig = ifelse(padj_male < 0.05, "Male padj < 0.05",
                        ifelse(padj_female < 0.05, "Female padj < 0.05",
                               "Neither significant")))

r <- cor(scatter_df$NES_male, scatter_df$NES_female, method = "spearman")

p_scatter <- ggplot(scatter_df, aes(x = NES_female, y = NES_male, color = sig)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_point(size = 2.6, alpha = 0.85) +
  scale_color_manual(values = c("Male padj < 0.05" = "#e31a1c",
                                "Female padj < 0.05" = "#1f78b4",
                                "Neither significant" = "grey55"),
                     name = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(title = "Hallmark NES - Male vs Female Input Cohorts",
       subtitle = sprintf("Spearman rho = %.2f across 50 sets", r),
       x = "Female NES", y = "Male NES")

ggsave(file.path(figures_dir, "13_gsea_hallmark_nes_scatter.png"),
       plot = p_scatter, width = 7.0, height = 6.0, dpi = 300, bg = "white")
message("  wrote 13_gsea_hallmark_nes_scatter.png")

# -----------------------------------------------------------------------------
# Fig E - NES distribution across all 50 sets, per contrast
# A well-behaved preranked run centres the NES distribution on 0. A median far
# below 0 means the ranking carries a global trend, not set-specific biology.
# -----------------------------------------------------------------------------
dist_df <- bind_rows(
  male_df   %>% mutate(contrast = "Male Input"),
  female_df %>% mutate(contrast = "Female Input")
)

med_tbl <- dist_df %>% group_by(contrast) %>%
  summarise(median_nes = median(NES))

p_dist <- ggplot(dist_df, aes(x = NES, fill = contrast)) +
  geom_histogram(binwidth = 0.2, colour = "grey40", alpha = 0.85) +
  geom_vline(data = med_tbl, aes(xintercept = median_nes, colour = contrast),
             linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values = c("Male Input" = "#e31a1c", "Female Input" = "#1f78b4"),
                    name = NULL) +
  scale_colour_manual(values = c("Male Input" = "#e31a1c", "Female Input" = "#1f78b4"),
                      guide = "none") +
  facet_grid(contrast ~ ., scales = "free_y") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)) +
  labs(title = "NES Distribution Across 50 Hallmark Sets",
       subtitle = sprintf("Dashed line = median NES. Male %.2f, Female %.2f. A null run centres on 0.",
                          med_tbl$median_nes[match("Male Input", med_tbl$contrast)],
                          med_tbl$median_nes[match("Female Input", med_tbl$contrast)]),
       x = "Normalized Enrichment Score (NES)", y = "Gene sets")

ggsave(file.path(figures_dir, "13_gsea_hallmark_nes_distribution.png"),
       plot = p_dist, width = 7.0, height = 6.0, dpi = 300, bg = "white")
message("  wrote 13_gsea_hallmark_nes_distribution.png")

message("\n13_gsea_figures.R completed. Figures in: ", figures_dir)
