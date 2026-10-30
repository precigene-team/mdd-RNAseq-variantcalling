# [[03_variant_calling/22_plot_burden_go_mf.R]]
# Goal: Control vs MDD burden for the Gene Ontology Molecular Function branch --
#       LEFT, the ten most impacted sets across all collections; RIGHT, the GO:MF
#       synaptic and neurotransmitter receptor sets. Tables from 16; styling 17.
#
# Reads:
#   16_burden_per_donor.csv
#   16_burden_summary.csv
#
# Outputs:
#   22_burden_go_mf.pdf / .png
#
# Usage:
#   Rscript 03_variant_calling/22_plot_burden_go_mf.R

suppressPackageStartupMessages(library(data.table))

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) curr <- dirname(curr)
  curr
}
source(file.path(get_project_root(), "03_variant_calling/17_plot_burden_control_vs_mdd.R"))

COLLECTION <- "GO_MF"
OUT_NAME   <- "22_burden_go_mf"

b <- load_burden()

pA <- top_panel(b$pd, b$sm, COLLECTION)
pB <- collection_panel(b$pd, b$sm, COLLECTION,
                       sprintf("B   %s synaptic / receptor function", COLLECTION),
                       "the receptor functions the brief asks about")

save_figure(compose_figure(pA, pB), b$ctx, OUT_NAME, MAX_PANEL_ROWS)



in_coll <- b$sm[collection == COLLECTION]
topA <- in_coll[genes_in_cohort >= MIN_GENES_PANEL][order(-abs(log2rr))][1:10]
message(sprintf("[22] %s: %d sets tested, panel A = ranks 1-10 by |effect|",
                COLLECTION, nrow(in_coll)))
for (i in seq_len(nrow(topA))) {
  message(sprintf("   A  %-44s ctrl=%4.0f mdd=%4.0f ratio=%.2f p=%.3f",
                  substr(pretty_name(topA$set_name[i]), 1, 44),
                  topA$median_ctrl[i], topA$median_mdd[i], topA$ratio[i], topA$p_glm[i]))
}
neuro <- in_coll[grepl(NEURO_PAT[[COLLECTION]], set_name)]
message(sprintf("[22] panel B: %d of %d sets match the synapse vocabulary",
                min(nrow(neuro), MAX_PANEL_ROWS), nrow(neuro)))
message("[22] done")
