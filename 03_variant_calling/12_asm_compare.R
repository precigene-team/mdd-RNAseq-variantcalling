# [[03_variant_calling/12_asm_compare.R]]
# Goal: Group-level comparison of allele-specific m6A (ASM) between MDD and Control.
#       Statistical analysis ONLY — the per-(locus,subject) allele-count matrix is built by
#       03_variant_calling/11_asm_build_matrix.py.
#
# Model (established machinery, per docs/asm_analysis_methodology.md §4):
#   Per locus: beta-binomial regression of IP allele counts on condition + sex, with
#   logit(Input allelic fraction) as offset (Input-as-null; 0.5 null when Input absent).
#   Primary engine: glmmTMB (betabinomial)  -> fallback lme4::glmer (binomial, (1|subject))
#   -> fallback stats::glm (quasibinomial). All are published, maintained implementations.
#   Cross-check: per-subject fisher.test(IP vs Input) log-OR, then wilcox.test across groups.
#   Multiple testing: p.adjust(..., "BH") per contrast. Sensitivity sweeps over skew/min-subjects.

suppressPackageStartupMessages({
  ok_glmmtmb <- requireNamespace("glmmTMB", quietly = TRUE)
  ok_lme4    <- requireNamespace("lme4", quietly = TRUE)
})

get_project_root <- function() {
  curr <- getwd()
  while (curr != "/" && !file.exists(file.path(curr, "run"))) {
    curr <- dirname(curr)
  }
  return(curr)
}
project_root <- get_project_root()

mat_file    <- file.path(project_root, "_data/RData/03_variant_calling/asm_matrix.tsv")
tables_dir  <- file.path(project_root, "results/03_variant_calling/tables")
fig_dir     <- file.path(project_root, "results/03_variant_calling/figures")
out_rds     <- file.path(project_root, "_data/RData/03_variant_calling/asm_results.rds")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Parameters ---
Q_CUT          <- 0.05
SKEW_CUT       <- 0.20
MIN_SUBJECTS   <- 8
MIN_IP_READS   <- 10
MIN_ALLELE     <- 2
MIN_INPUT_READS<- 8
FDR_CUT        <- 0.05
SWEEP_SKEW     <- c(0.1, 0.2, 0.3)
SWEEP_MINSUBJ  <- c(6, 8, 12)

logit <- function(p) log(p / (1 - p))

# --- Load matrix ---
message("[INFO] Loading ASM matrix: ", mat_file)
mat <- read.delim(mat_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))
if (!nrow(mat)) stop("ASM matrix is empty; run 11_asm_build_matrix.py first.")
input_available <- any(mat$in_avail == 1)
message("[INFO] rows: ", nrow(mat), "; input counts available: ", input_available)

mat$condition <- factor(mat$condition, levels = c("ctrl", "mdd"))
mat$sex       <- factor(mat$sex, levels = c("male", "female"))
mat$ip_b      <- mat$ip_total - mat$ip_ref
mat$ip_af     <- mat$ip_ref / mat$ip_total
mat$skew      <- abs(mat$ip_af - 0.5)
# Input fraction / offset (smoothed). 0.5-null when Input counts absent.
mat$in_af <- ifelse(mat$in_avail == 1,
                    (mat$in_ref + 0.5) / (mat$in_total + 1),
                    0.5)
mat$off   <- logit(pmax(pmin(mat$in_af, 1 - 1e-6), 1e-6))
mat$has_input <- mat$in_avail == 1 & mat$in_total >= MIN_INPUT_READS

# --- Observation inclusion (for a given skew cutoff) ---
include_obs <- function(d, skew_cut) {
  keep <- d$ip_total >= MIN_IP_READS &
          d$ip_ref >= MIN_ALLELE &
          d$ip_b >= MIN_ALLELE &
          d$min_q < Q_CUT &
          d$skew >= skew_cut
  if (input_available) keep <- keep & d$has_input
  keep
}

# --- Model fitter with engine fallback. Returns list(coef, se, p, phi) for a term. ---
# adj_sex = TRUE adds sex as a covariate to the "condition" formula (pooled primary);
# within-sex sub-contrasts must keep adj_sex = FALSE (sex is constant in the subset).
fit_locus <- function(d, term, offset_ok = TRUE, adj_sex = FALSE) {
  eff <- if (adj_sex && term == "condition") "condition + sex" else term
  formula_glmmtmb <- as.formula(sprintf("cbind(ip_ref, ip_b) ~ %s + offset(off)", eff))
  formula_glmer   <- as.formula(sprintf("cbind(ip_ref, ip_b) ~ %s + (1|subject) + offset(off)", eff))
  formula_glm     <- as.formula(sprintf("cbind(ip_ref, ip_b) ~ %s + offset(off)", eff))
  coef_name <- switch(term,
                      "condition" = "conditionmdd",
                      "condition * sex" = "conditionmdd:sexfemale",
                      "condition*sex"   = "conditionmdd:sexfemale")
  if (is.null(coef_name)) stop("unhandled term: ", term)
  fit <- NULL
  engine <- NULL
  if (ok_glmmtmb) {
    fit <- tryCatch(glmmTMB::glmmTMB(formula_glmmtmb, family = glmmTMB::betabinomial(),
                                     data = d), error = function(e) NULL)
    if (!is.null(fit)) engine <- "glmmTMB"
  }
  if (is.null(fit) && ok_lme4) {
    fit <- tryCatch(lme4::glmer(formula_glmer, family = binomial, data = d),
                    error = function(e) NULL)
    if (!is.null(fit)) engine <- "lme4"
  }
  if (is.null(fit)) {
    fit <- tryCatch(glm(formula_glm, family = quasibinomial, data = d),
                    error = function(e) NULL)
    if (!is.null(fit)) engine <- "glm"
  }
  if (is.null(fit)) return(list(coef = NA, se = NA, p = NA, phi = NA, engine = NA))
  s <- summary(fit)
  if (engine == "glmmTMB") {
    ct <- s$coefficients$cond
  } else {
    ct <- s$coefficients
  }
  if (!coef_name %in% rownames(ct)) {
    # fall back to the last row of the fixed-effect table
    coef_name <- tail(rownames(ct), 1)
  }
  phi <- NA_real_
  if (engine == "glm") phi <- s$dispersion
  list(coef = ct[coef_name, 1], se = ct[coef_name, 2],
       p = ct[coef_name, ncol(ct)], phi = phi, engine = engine)
}

# --- Per-locus summary for a given (skew_cut, min_subjects) ---
analyze <- function(skew_cut, min_subjects) {
  d <- mat[include_obs(mat, skew_cut), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  by_locus <- split(d, d$locus_id)
  n_snps <- sapply(by_locus, nrow)
  loci <- names(by_locus)[n_snps >= min_subjects]
  if (!length(loci)) return(NULL)
  res <- list()
  for (id in loci) {
    dd <- by_locus[[id]]
    n_mdd <- sum(dd$condition == "mdd")
    n_ctrl <- sum(dd$condition == "ctrl")
    if (n_mdd < 2 || n_ctrl < 2) next
    prim <- fit_locus(dd, "condition", adj_sex = length(unique(as.character(dd$sex))) > 1)
    # within-sex (exploratory)
    dm <- dd[dd$sex == "male", , drop = FALSE]
    df_ <- dd[dd$sex == "female", , drop = FALSE]
    male <- if (sum(dm$condition == "mdd") >= 2 && sum(dm$condition == "ctrl") >= 2)
      fit_locus(dm, "condition") else list(coef = NA, se = NA, p = NA, phi = NA, engine = NA)
    female <- if (sum(df_$condition == "mdd") >= 2 && sum(df_$condition == "ctrl") >= 2)
      fit_locus(df_, "condition") else list(coef = NA, se = NA, p = NA, phi = NA, engine = NA)
    inter <- fit_locus(dd, "condition * sex")
    # cross-check: per-subject Fisher log-OR (IP vs Input), then Wilcoxon between groups
    mwu_p <- NA_real_
    fisher_n_sig <- NA_integer_
    if (input_available && all(dd$has_input)) {
      lor <- rep(NA_real_, nrow(dd))
      fsig <- 0L
      for (i in seq_len(nrow(dd))) {
        tb <- matrix(c(dd$ip_ref[i], dd$ip_b[i], dd$in_ref[i], dd$in_total[i] - dd$in_ref[i]),
                     nrow = 2, byrow = TRUE)
        ft <- tryCatch(fisher.test(tb), error = function(e) NULL)
        if (is.null(ft) || is.null(ft$p.value)) next
        if (ft$p.value < 0.05) fsig <- fsig + 1L
        est <- unname(ft$estimate)
        if (is.finite(est) && est > 0) lor[i] <- log(est)
      }
      fisher_n_sig <- fsig
      g <- dd$condition
      keep <- is.finite(lor)
      if (sum(keep & g == "mdd") >= 3 && sum(keep & g == "ctrl") >= 3) {
        mwu_p <- tryCatch(wilcox.test(lor[keep & g == "mdd"], lor[keep & g == "ctrl"])$p.value,
                          error = function(e) NA_real_)
      }
    }
    res[[id]] <- data.frame(
      locus_id = id, chr = dd$chr[1], start = dd$start[1], end = dd$end[1],
      gene = dd$gene[1], n_subjects = nrow(dd), n_mdd = n_mdd, n_ctrl = n_ctrl,
      n_snps_total = sum(dd$n_snps), ip_total = sum(dd$ip_total),
      mdd_af = mean(dd$ip_af[dd$condition == "mdd"]), ctrl_af = mean(dd$ip_af[dd$condition == "ctrl"]),
      coef = prim$coef, se = prim$se, p = prim$p, engine = prim$engine,
      male_coef = male$coef, male_p = male$p, male_engine = male$engine,
      female_coef = female$coef, female_p = female$p, female_engine = female$engine,
      interaction_coef = inter$coef, interaction_p = inter$p, interaction_engine = inter$engine,
      mwu_p = mwu_p, fisher_n_sig = fisher_n_sig,
      stringsAsFactors = FALSE
    )
  }
  if (!length(res)) return(NULL)
  out <- do.call(rbind, res)
  out$fdr <- p.adjust(out$p, method = "BH")
  out$male_fdr <- p.adjust(out$male_p, method = "BH")
  out$female_fdr <- p.adjust(out$female_p, method = "BH")
  if (input_available) out$mwu_fdr <- p.adjust(out$mwu_p, method = "BH")
  # Confidence grade for the primary contrast: High requires balanced groups AND a
  # sex-specific FDR<0.05 signal; otherwise Low (control-poor / no sex-replicated hit).
  sex_hit <- !is.na(out$male_fdr) & out$male_fdr < FDR_CUT |
             !is.na(out$female_fdr) & out$female_fdr < FDR_CUT
  out$confidence <- ifelse(out$n_ctrl >= 5 & sex_hit, "High", "Low")
  out
}

# --- Main run (default thresholds) ---
message("[INFO] Primary analysis: skew=", SKEW_CUT, ", min_subjects=", MIN_SUBJECTS)
res <- analyze(SKEW_CUT, MIN_SUBJECTS)
if (is.null(res)) stop("No loci passed filters at the default thresholds.")
message("[INFO] loci tested: ", nrow(res), "; FDR<", FDR_CUT, ": ", sum(res$fdr < FDR_CUT))

locus_table <- file.path(tables_dir, "12_asm_locus_table.csv")
sig_table   <- file.path(tables_dir, "12_asm_primary_significant.csv")
write.csv(res, locus_table, row.names = FALSE)
write.csv(res[res$fdr < FDR_CUT & !is.na(res$fdr), , drop = FALSE], sig_table, row.names = FALSE)
message("[INFO] wrote: ", locus_table, " / ", sig_table)

if (input_available) {
  cross <- res[!is.na(res$mwu_p), c("locus_id", "gene", "coef", "p", "fdr", "mwu_p", "mwu_fdr", "fisher_n_sig")]
  cross$rank_agree <- cor(res$p, res$mwu_p, method = "spearman", use = "complete.obs")
  write.csv(cross, file.path(tables_dir, "12_asm_crosscheck.csv"), row.names = FALSE)
  message("[INFO] Spearman rank corr (GLM p vs MWU p): ", round(cross$rank_agree[1], 3))
}

# --- Sensitivity sweep ---
sweep <- expand.grid(skew_cut = SWEEP_SKEW, min_subjects = SWEEP_MINSUBJ, KEEP.OUT.ATTRS = FALSE)
sweep$n_loci <- sweep$n_primary_sig <- sweep$n_male_sig <- sweep$n_female_sig <- NA_integer_
for (i in seq_len(nrow(sweep))) {
  r <- analyze(sweep$skew_cut[i], sweep$min_subjects[i])
  if (!is.null(r)) {
    sweep$n_loci[i] <- nrow(r)
    sweep$n_primary_sig[i] <- sum(r$fdr < FDR_CUT, na.rm = TRUE)
    sweep$n_male_sig[i] <- sum(r$male_fdr < FDR_CUT, na.rm = TRUE)
    sweep$n_female_sig[i] <- sum(r$female_fdr < FDR_CUT, na.rm = TRUE)
  }
}
write.csv(sweep, file.path(tables_dir, "12_asm_sensitivity.csv"), row.names = FALSE)
message("[INFO] wrote sensitivity: ", file.path(tables_dir, "12_asm_sensitivity.csv"))

# --- Figures (base R) ---
pdf(file.path(fig_dir, "12_asm_volcano.pdf"), width = 6, height = 6)
x <- res$coef
y <- -log10(res$fdr)
plot(x, y, pch = 20, col = ifelse(res$fdr < FDR_CUT, "firebrick", "grey60"),
     xlab = "condition log-OR (IP allelic fraction, MDD vs Ctrl)",
     ylab = "-log10(FDR)", main = paste0("ASM group comparison (n loci = ", nrow(res), ")"))
abline(h = -log10(FDR_CUT), lty = 2, col = "grey40")
legend("topright", legend = c("FDR<0.05", "ns"), pch = 20, col = c("firebrick", "grey60"))
dev.off()

pdf(file.path(fig_dir, "12_asm_sensitivity.pdf"), width = 7, height = 5)
bp <- matrix(sweep$n_primary_sig, nrow = length(SWEEP_SKEW), byrow = TRUE,
             dimnames = list(paste0("skew=", SWEEP_SKEW), paste0("minSubj=", SWEEP_MINSUBJ)))
barplot(bp, beside = TRUE, legend.text = TRUE, col = c("steelblue", "grey70", "firebrick"),
        xlab = "minimum subjects per locus", ylab = "significant loci (FDR<0.05)",
        main = "ASM sensitivity: significant loci")
dev.off()
message("[INFO] wrote figures to ", fig_dir)

saveRDS(list(results = res, sensitivity = sweep, params = list(
  q_cut = Q_CUT, skew_cut = SKEW_CUT, min_subjects = MIN_SUBJECTS, min_ip_reads = MIN_IP_READS,
  min_allele = MIN_ALLELE, min_input_reads = MIN_INPUT_READS, fdr_cut = FDR_CUT,
  input_available = input_available)), out_rds)
message("[INFO] wrote RDS: ", out_rds)
message("[INFO] done.")
