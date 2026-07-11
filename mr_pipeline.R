# ============================================================
# Mendelian randomisation pipeline: 14-protein panel <-> lung
# cancer <-> smoking behaviour. Run sections in order; each
# writes files the next section reads. Expects a working
# directory containing data/ (raw inputs) and results/ (outputs).
# ============================================================

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("MRCIEU/TwoSampleMR")
}
if (!requireNamespace("MendelianRandomization", quietly = TRUE)) install.packages("MendelianRandomization")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
library(data.table)
library(ggplot2)
library(patchwork)

DATA_DIR <- "data"
OUT_DIR  <- "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

GENES <- c("ALPP","CDCP1","CEACAM5","CXCL17","GDF15","LAMP3","MMP12",
           "PIGR","PLAUR","PRSS8","SFTPA1","SFTPD","TNFSF13B","WFDC2")

UKBPPP_N <- 34557
SMKINIT_NCASE <- 140856     # ever-smokers, GSCAN excl. UKB/23andMe
SMKINIT_NCONTROL <- 108315  # never-smokers

DECOMPRESS <- if (nzchar(Sys.which("pigz"))) "pigz -dc" else "gzip -dc"

# ---- Shared helpers, used across multiple sections below ----

get_slot_safe <- function(obj, ...) {
  if (is.null(obj)) return(NA)
  avail <- methods::slotNames(obj)
  for (cand in c(...)) if (cand %in% avail) return(methods::slot(obj, cand))
  NA
}

# Greedy LD clumping (r2<=0.001): sorts by significance, keeps top SNP as
# independent, removes correlated SNPs, repeats.
greedy_clump <- function(rsid, sig, ld_mat, ld_labels, r2_threshold = 0.001, sig_is_log10 = FALSE) {
  keep <- rsid %in% ld_labels
  if (sum(keep) == 0) return(character(0))
  rsid <- rsid[keep]; sig <- sig[keep]
  idx_in_mat <- match(rsid, ld_labels)
  order_by_sig <- if (sig_is_log10) order(-sig) else order(sig)
  independent <- character(0)
  eligible <- rep(TRUE, length(rsid))
  for (i in order_by_sig) {
    if (!eligible[i]) next
    independent <- c(independent, rsid[i])
    r2_to_i <- ld_mat[idx_in_mat[i], idx_in_mat]^2
    eligible[which(r2_to_i > r2_threshold)] <- FALSE
    eligible[i] <- FALSE
  }
  independent
}

# Reads an LD matrix, detecting whether SNP labels are the first column or
# the column names themselves.
read_ld_matrix <- function(ld_file) {
  ld_raw <- fread(ld_file, header = TRUE)
  first_col_is_labels <- !is.numeric(ld_raw[[1]]) ||
    !all(colnames(ld_raw) %in% as.character(seq_len(ncol(ld_raw))))
  if (first_col_is_labels && colnames(ld_raw)[1] %in% c("V1","SNP","rsid","","id")) {
    ld_labels <- ld_raw[[1]]
    ld_mat <- as.matrix(ld_raw[, -1, with = FALSE])
  } else {
    ld_labels <- colnames(ld_raw)
    ld_mat <- as.matrix(ld_raw)
  }
  rownames(ld_mat) <- ld_labels
  list(mat = ld_mat, labels = ld_labels)
}

radial_outlier_test <- function(rsid, bx, by, byse) {
  w <- bx^2 / byse^2
  ratio <- by / bx
  beta_ivw <- sum(w * ratio) / sum(w)
  q_i <- w * (ratio - beta_ivw)^2
  p_i <- pchisq(q_i, df = 1, lower.tail = FALSE)
  data.table(rsid = rsid, q_contribution = q_i, p_outlier = p_i,
             is_outlier = p_i < (0.05 / length(rsid)))
}

run_ivw_manual <- function(bx, by, byse) {
  w <- bx^2 / byse^2
  beta <- sum(bx * by / byse^2) / sum(w)
  se <- sqrt(1 / sum(w))
  list(beta = beta, se = se, p = 2*pnorm(-abs(beta/se)), n_snps = length(bx))
}

fmt_p <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.3f", p)))

# ============================================================
# SECTION 1: extract cis-pQTL candidate instruments (14 proteins)
# ============================================================

WINDOW_KB <- 1000; GWAS_SIG_P <- 5e-8; INFO_MIN <- 0.8; MAF_MIN <- 0.01

gene_coords <- data.table(
  gene  = GENES,
  chrom = c("2","3","19","19","19","3","11","1","19","16","10","10","13","20"),
  tss   = c(232382889, 45146483, 41730519, 42443029, 18389176, 183163839,
            102874990, 206949367, 43670547, 31135727, 79649708, 79983996,
            108308484, 45481533)
)

pqtl_tar_files <- list.files(file.path(DATA_DIR, "pqtl_raw"), pattern = "\\.tar$", full.names = TRUE)
names(pqtl_tar_files) <- toupper(sub("_.*", "", basename(pqtl_tar_files)))

snp_map_cache <- new.env()
get_snp_map <- function(chrom) {
  key <- as.character(chrom)
  if (!is.null(snp_map_cache[[key]])) return(snp_map_cache[[key]])
  map_file <- file.path(DATA_DIR, "snp_map", paste0("olink_rsid_map_mac5_info03_b0_7_chr", chrom, "_patched_v2.tsv.gz"))
  if (!file.exists(map_file)) return(NULL)
  map <- fread(map_file)
  snp_map_cache[[key]] <- map
  map
}

extract_cis_window <- function(gene) {
  if (!gene %in% names(pqtl_tar_files)) return(NULL)
  g_row <- gene_coords[gene_coords$gene == gene, ]
  chrom <- g_row$chrom; tss <- g_row$tss

  contents <- untar(pqtl_tar_files[gene], list = TRUE)
  target <- grep(paste0("discovery_chr", chrom, "_"), contents, value = TRUE)
  if (length(target) == 0) return(NULL)
  tmp_dir <- tempfile(); dir.create(tmp_dir)
  untar(pqtl_tar_files[gene], files = target, exdir = tmp_dir)
  pqtl <- fread(file.path(tmp_dir, target))
  pqtl[, PVAL := 10^(-LOG10P)]

  snp_map <- get_snp_map(chrom)
  if (is.null(snp_map)) return(NULL)
  merged <- merge(pqtl, snp_map, by = "ID", all.x = TRUE)

  window <- merged[POS38 >= (tss - WINDOW_KB*1000) & POS38 <= (tss + WINDOW_KB*1000)]
  window[, qc_pass := INFO >= INFO_MIN & pmin(A1FREQ, 1 - A1FREQ) >= MAF_MIN]

  fwrite(window, file.path(DATA_DIR, "cis_windows", paste0(gene, "_cis_window_full.csv.gz")))
  fwrite(window[qc_pass == TRUE & PVAL < GWAS_SIG_P],
         file.path(DATA_DIR, "cis_windows", paste0(gene, "_instrument_candidates.csv")))
  data.frame(gene = gene, n_total = nrow(window), n_qc_pass = sum(window$qc_pass))
}

dir.create(file.path(DATA_DIR, "cis_windows"), showWarnings = FALSE, recursive = TRUE)
cis_summary <- do.call(rbind, lapply(GENES, extract_cis_window))
print(cis_summary)

# ============================================================
# SECTION 2: clump protein instruments (r2<=0.001, F>=10)
# ============================================================

F_STAT_MIN <- 10  # CHISQ = F-statistic for a single-variant instrument

clump_protein <- function(gene) {
  cand_file <- file.path(DATA_DIR, "cis_windows", paste0(gene, "_instrument_candidates.csv"))
  ld_file <- file.path(DATA_DIR, "ld", "protein", paste0(gene, "_ld.csv"))
  if (!file.exists(cand_file) || !file.exists(ld_file)) return(NULL)

  cand <- fread(cand_file)
  cand <- cand[!grepl("^[0-9XY]+:[0-9]+_", rsid) & !is.na(rsid)]
  if ("CHISQ" %in% colnames(cand)) cand <- cand[CHISQ >= F_STAT_MIN]
  if (nrow(cand) == 0) return(data.frame(gene = gene, n_clumped = 0))

  ld <- read_ld_matrix(ld_file)
  independent_rsids <- greedy_clump(cand$rsid, cand$LOG10P, ld$mat, ld$labels, sig_is_log10 = TRUE)
  final <- cand[rsid %in% independent_rsids]
  fwrite(final, file.path(DATA_DIR, "instruments", "protein", paste0(gene, "_final_instruments.csv")))
  data.frame(gene = gene, n_candidates = nrow(cand), n_clumped = length(independent_rsids))
}

dir.create(file.path(DATA_DIR, "instruments", "protein"), showWarnings = FALSE, recursive = TRUE)
clump_summary <- do.call(rbind, lapply(GENES, clump_protein))
print(clump_summary)

# ============================================================
# SECTION 3: extract lung cancer GWAS rows matching protein instruments
# ============================================================

LC_OUTCOMES <- c("lung_cancer", "lung_adenocarcinoma", "lung_squamous",
                  "lung_eversmoker", "lung_neversmoker")

all_protein_rsids <- unique(unlist(lapply(GENES, function(g) {
  f <- file.path(DATA_DIR, "instruments", "protein", paste0(g, "_final_instruments.csv"))
  if (file.exists(f)) fread(f)$rsid else character(0)
})))

rsid_file <- tempfile(fileext = ".txt")
writeLines(all_protein_rsids, rsid_file)
awk_lc <- tempfile(fileext = ".awk")
writeLines(c('NR==FNR { rsid[$1]=1; next }', 'FNR==1 { next }',
             '($2 in rsid) || ($13 in rsid) { print }'), awk_lc)

dir.create(file.path(DATA_DIR, "lung_gwas", "extracted"), showWarnings = FALSE, recursive = TRUE)
for (outcome in LC_OUTCOMES) {
  h_file <- list.files(file.path(DATA_DIR, "lung_gwas", outcome, "harmonized"),
                        pattern = "\\.h\\.tsv\\.gz$", full.names = TRUE)[1]
  if (is.na(h_file)) next
  header <- readLines(h_file, n = 1)
  out_file <- file.path(DATA_DIR, "lung_gwas", "extracted", paste0(outcome, "_matched.tsv"))
  writeLines(header, out_file)
  cmd <- paste(DECOMPRESS, shQuote(h_file), "| awk -f", shQuote(awk_lc), shQuote(rsid_file), "-", ">>", shQuote(out_file))
  system(cmd)
}

# ============================================================
# SECTION 4: MR, protein -> lung cancer (IVW/Wald, Egger, WME, MBE)
# ============================================================

lc_outcome_meta <- data.frame(
  outcome  = LC_OUTCOMES,
  ncase    = c(29266, 11273, 7426, 23223, 2355),
  ncontrol = c(56450, 55483, 55627, 16964, 7504)
)
lc_outcome_meta$samplesize <- lc_outcome_meta$ncase + lc_outcome_meta$ncontrol

mr_protein_to_lc <- function(gene, outcome_row) {
  outcome <- outcome_row$outcome
  cand_file <- file.path(DATA_DIR, "instruments", "protein", paste0(gene, "_final_instruments.csv"))
  out_file  <- file.path(DATA_DIR, "lung_gwas", "extracted", paste0(outcome, "_matched.tsv"))
  if (!file.exists(cand_file) || !file.exists(out_file)) return(data.frame(gene = gene, outcome = outcome, status = "missing input"))

  exp_raw <- fread(cand_file); out_raw <- fread(out_file)
  if (nrow(exp_raw) == 0) return(data.frame(gene = gene, outcome = outcome, status = "no candidates"))

  exp_dat <- TwoSampleMR::format_data(as.data.frame(exp_raw), type = "exposure",
    snp_col = "rsid", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "PVAL")
  exp_dat$samplesize.exposure <- UKBPPP_N; exp_dat$exposure <- gene

  out_df <- as.data.frame(out_raw)
  out_df$ncase <- outcome_row$ncase; out_df$ncontrol <- outcome_row$ncontrol; out_df$samplesize <- outcome_row$samplesize
  out_dat <- TwoSampleMR::format_data(out_df, type = "outcome",
    snp_col = "hm_rsid", beta_col = "hm_beta", se_col = "standard_error",
    effect_allele_col = "hm_effect_allele", other_allele_col = "hm_other_allele",
    eaf_col = "hm_effect_allele_frequency", pval_col = "p_value",
    ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize")
  out_dat$outcome <- outcome

  harm <- tryCatch(TwoSampleMR::harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
  if (is.null(harm) || nrow(harm) == 0) return(data.frame(gene = gene, outcome = outcome, status = "harmonisation failed"))
  harm <- harm[harm$mr_keep, ]
  if (nrow(harm) == 0) return(data.frame(gene = gene, outcome = outcome, status = "all SNPs dropped"))

  harm <- tryCatch(TwoSampleMR::steiger_filtering(harm), error = function(e) { harm$steiger_dir <- NA; harm })
  harm_steiger <- harm[!is.na(harm$steiger_dir) & harm$steiger_dir, ]
  if (nrow(harm_steiger) == 0) return(data.frame(gene = gene, outcome = outcome, status = "failed Steiger"))

  mrin <- MendelianRandomization::mr_input(bx = harm_steiger$beta.exposure, bxse = harm_steiger$se.exposure,
                                            by = harm_steiger$beta.outcome, byse = harm_steiger$se.outcome)
  n_snps <- nrow(harm_steiger)
  method <- if (n_snps == 1) "wald_ratio" else "ivw"
  ivw <- MendelianRandomization::mr_ivw(mrin)
  egger <- NULL; wme <- NULL; mbe <- NULL
  if (n_snps >= 3) {
    egger <- tryCatch(MendelianRandomization::mr_egger(mrin), error = function(e) NULL)
    wme   <- tryCatch(MendelianRandomization::mr_median(mrin, weighting = "weighted"), error = function(e) NULL)
    mbe   <- tryCatch(MendelianRandomization::mr_mbe(mrin), error = function(e) NULL)
  }
  het_p <- if (n_snps == 1) NA else tryCatch(ivw@Heter.Stat[2], error = function(e) NA)

  data.frame(
    gene = gene, outcome = outcome, status = "ok", method = method, n_snps = n_snps,
    beta = get_slot_safe(ivw, "Estimate"), se = get_slot_safe(ivw, "StdError"), p = get_slot_safe(ivw, "Pvalue"),
    OR = exp(get_slot_safe(ivw, "Estimate")), OR_lo = exp(get_slot_safe(ivw, "CILower")), OR_hi = exp(get_slot_safe(ivw, "CIUpper")),
    heterogeneity_p = het_p, egger_intercept = get_slot_safe(egger, "Intercept"), egger_intercept_p = get_slot_safe(egger, "Pvalue.Int"),
    wme_beta = get_slot_safe(wme, "Estimate"), wme_se = get_slot_safe(wme, "StdError"), wme_p = get_slot_safe(wme, "Pvalue"),
    mbe_beta = get_slot_safe(mbe, "Estimate"), mbe_se = get_slot_safe(mbe, "StdError"), mbe_p = get_slot_safe(mbe, "Pvalue")
  )
}

lc_results <- list()
for (g in GENES) for (i in seq_len(nrow(lc_outcome_meta))) lc_results[[paste(g, i)]] <- mr_protein_to_lc(g, lc_outcome_meta[i, ])
lc_results_df <- rbindlist(lc_results, fill = TRUE)
fwrite(lc_results_df, file.path(OUT_DIR, "mr_protein_to_lung_cancer.csv"))

# ---- Table 1 ----
lc_valid <- lc_results_df[status == "ok"]
lc_valid[, P_FDR := p.adjust(p, method = "BH")]
outcome_labels <- c(lung_cancer = "Lung cancer (overall)", lung_adenocarcinoma = "Lung adenocarcinoma",
                     lung_squamous = "Lung squamous cell carcinoma", lung_eversmoker = "Lung cancer, ever-smokers",
                     lung_neversmoker = "Lung cancer, never-smokers")
lc_valid[, Outcome := outcome_labels[outcome]]
lc_valid[, Method := ifelse(method == "wald_ratio", "Wald ratio", "IVW")]
lc_valid[, `OR (95% CI)` := sprintf("%.2f (%.2f-%.2f)", OR, OR_lo, OR_hi)]
lc_valid[, `WME (95% CI)` := ifelse(is.na(wme_beta), "-", sprintf("%.3f (%.3f-%.3f)", wme_beta, wme_beta-1.96*wme_se, wme_beta+1.96*wme_se))]
lc_valid[, `MBE (95% CI)` := ifelse(is.na(mbe_beta), "-", sprintf("%.3f (%.3f-%.3f)", mbe_beta, mbe_beta-1.96*mbe_se, mbe_beta+1.96*mbe_se))]

table1 <- lc_valid[, .(
  Protein = gene, Outcome, Method, nSNP = n_snps, `OR (95% CI)`,
  P = fmt_p(p), `P (FDR-adjusted)` = fmt_p(P_FDR),
  `Weighted median (beta, 95% CI)` = `WME (95% CI)`, `WME P` = fmt_p(wme_p),
  `Weighted mode (beta, 95% CI)` = `MBE (95% CI)`, `MBE P` = fmt_p(mbe_p),
  `Heterogeneity P` = ifelse(is.na(heterogeneity_p), "-", sprintf("%.3f", heterogeneity_p)),
  `Egger intercept P` = ifelse(is.na(egger_intercept_p), "-", sprintf("%.3f", egger_intercept_p))
)]
outcome_order <- c("Lung cancer (overall)", "Lung adenocarcinoma", "Lung squamous cell carcinoma",
                    "Lung cancer, ever-smokers", "Lung cancer, never-smokers")
table1[, Outcome := factor(Outcome, levels = outcome_order)]
setorder(table1, Protein, Outcome)
table1[, Outcome := as.character(Outcome)]

con <- file(file.path(OUT_DIR, "Supplementary_Table_1.tsv"), "w")
writeLines(c("Supplementary Table 1. MR estimates for 14 panel proteins and lung cancer risk.",
             "Cis-pQTLs from UKB-PPP (n=34,557; P<5e-8, MAF>0.01, INFO>0.8, F>10, LD r2<=0.001) vs ILCCO/TRICL lung cancer (UKB excluded). IVW/Wald, weighted median/mode, MR-Egger for >=3 SNPs. FDR (BH) across all tests.", ""), con)
close(con)
fwrite(table1, file.path(OUT_DIR, "Supplementary_Table_1.tsv"), append = TRUE, col.names = TRUE, sep = "\t")

# ============================================================
# SECTION 5: MR, protein -> smoking (collider-bias/pleiotropy check)
# ============================================================

smoking_files <- c(SmokingInitiation = "SmokingInitiation.WithoutUKB.txt.gz",
                    CigarettesPerDay  = "CigarettesPerDay.WithoutUKB.txt.gz")

protein_instruments_all <- rbindlist(lapply(GENES, function(g) {
  f <- file.path(DATA_DIR, "instruments", "protein", paste0(g, "_final_instruments.csv"))
  if (!file.exists(f)) return(NULL)
  d <- fread(f); d$gene <- g; d
}), fill = TRUE)

rsid_file2 <- tempfile(fileext = ".txt")
writeLines(unique(protein_instruments_all$rsid), rsid_file2)
smoking_header <- strsplit(system(paste(DECOMPRESS, shQuote(file.path(DATA_DIR, "smoking_gwas", smoking_files[1])), "| head -1"), intern = TRUE), "\t")[[1]]
rsid_col_idx <- which(toupper(smoking_header) == "RSID")
awk_smk <- tempfile(fileext = ".awk")
writeLines(c('NR==FNR { rsid[$1]=1; next }', 'FNR==1 { next }',
             sprintf('($%d in rsid) { print }', rsid_col_idx)), awk_smk)

smoking_dat <- list()
for (pheno in names(smoking_files)) {
  cmd <- paste(DECOMPRESS, shQuote(file.path(DATA_DIR, "smoking_gwas", smoking_files[pheno])), "| awk -f", shQuote(awk_smk), shQuote(rsid_file2), "-")
  d <- fread(text = paste(system(cmd, intern = TRUE), collapse = "\n"), header = FALSE, sep = "\t")
  setnames(d, seq_along(smoking_header), smoking_header)
  d$phenotype <- pheno
  smoking_dat[[pheno]] <- d
}
smoking_dat <- rbindlist(smoking_dat, fill = TRUE)

mr_protein_to_smoking <- function(target_gene, target_phenotype) {
  exp_raw <- protein_instruments_all[gene == target_gene]
  out_raw <- smoking_dat[phenotype == target_phenotype]
  if (nrow(exp_raw) == 0 || nrow(out_raw) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "no data"))

  exp_dat <- TwoSampleMR::format_data(as.data.frame(exp_raw), type = "exposure",
    snp_col = "rsid", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf_col = "A1FREQ", pval_col = "PVAL")
  exp_dat$samplesize.exposure <- UKBPPP_N; exp_dat$exposure <- target_gene

  out_df <- as.data.frame(out_raw)
  is_binary <- target_phenotype == "SmokingInitiation"
  if (is_binary) {
    out_df$ncase <- SMKINIT_NCASE; out_df$ncontrol <- SMKINIT_NCONTROL; out_df$samplesize <- SMKINIT_NCASE + SMKINIT_NCONTROL
    out_dat <- TwoSampleMR::format_data(out_df, type = "outcome",
      snp_col = "RSID", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF",
      eaf_col = "AF", pval_col = "PVALUE", ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize")
  } else {
    out_df$samplesize <- out_df$N
    out_dat <- TwoSampleMR::format_data(out_df, type = "outcome",
      snp_col = "RSID", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF",
      eaf_col = "AF", pval_col = "PVALUE", samplesize_col = "samplesize")
  }
  out_dat$outcome <- target_phenotype

  harm <- tryCatch(TwoSampleMR::harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
  if (is.null(harm) || nrow(harm) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "harmonisation failed"))
  harm <- harm[harm$mr_keep, ]
  if (nrow(harm) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "all SNPs dropped"))

  harm_steiger <- tryCatch(TwoSampleMR::steiger_filtering(harm), error = function(e) NULL)
  if (is.null(harm_steiger)) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "steiger failed"))
  harm_final <- harm_steiger[!is.na(harm_steiger$steiger_dir) & harm_steiger$steiger_dir, ]
  if (nrow(harm_final) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "failed Steiger"))

  n_snps <- nrow(harm_final)
  mrin <- MendelianRandomization::mr_input(bx = harm_final$beta.exposure, bxse = harm_final$se.exposure,
                                            by = harm_final$beta.outcome, byse = harm_final$se.outcome)
  ivw <- MendelianRandomization::mr_ivw(mrin)
  egger <- NULL; wme <- NULL; mbe <- NULL
  if (n_snps >= 3) {
    egger <- tryCatch(MendelianRandomization::mr_egger(mrin), error = function(e) NULL)
    wme   <- tryCatch(MendelianRandomization::mr_median(mrin, weighting = "weighted"), error = function(e) NULL)
    mbe   <- tryCatch(MendelianRandomization::mr_mbe(mrin), error = function(e) NULL)
  }
  het_p <- if (n_snps == 1) NA else tryCatch(ivw@Heter.Stat[2], error = function(e) NA)
  ivw_beta <- get_slot_safe(ivw, "Estimate")

  data.table(
    gene = target_gene, phenotype = target_phenotype, status = "ok",
    outcome_type = ifelse(is_binary, "binary", "continuous"), n_snps = n_snps,
    beta = ivw_beta, se = get_slot_safe(ivw, "StdError"), p = get_slot_safe(ivw, "Pvalue"),
    OR_or_NA = ifelse(is_binary, exp(ivw_beta), NA),
    OR_lo_or_NA = ifelse(is_binary, exp(get_slot_safe(ivw, "CILower")), NA),
    OR_hi_or_NA = ifelse(is_binary, exp(get_slot_safe(ivw, "CIUpper")), NA),
    beta_lo_if_continuous = ifelse(!is_binary, get_slot_safe(ivw, "CILower"), NA),
    beta_hi_if_continuous = ifelse(!is_binary, get_slot_safe(ivw, "CIUpper"), NA),
    heterogeneity_p = het_p, egger_intercept_p = get_slot_safe(egger, "Pvalue.Int"),
    wme_beta = get_slot_safe(wme, "Estimate"), wme_se = get_slot_safe(wme, "StdError"), wme_p = get_slot_safe(wme, "Pvalue"),
    mbe_beta = get_slot_safe(mbe, "Estimate"), mbe_se = get_slot_safe(mbe, "StdError"), mbe_p = get_slot_safe(mbe, "Pvalue")
  )
}

p2s_phenotypes <- c("SmokingInitiation", "CigarettesPerDay")
p2s_results <- rbindlist(lapply(GENES, function(g) rbindlist(lapply(p2s_phenotypes, function(p) mr_protein_to_smoking(g, p)), fill = TRUE)), fill = TRUE)
fwrite(p2s_results, file.path(OUT_DIR, "mr_protein_to_smoking.csv"))

# ---- Table 2 (FDR across all tests combined) ----
p2s_valid <- p2s_results[status == "ok"]
p2s_valid[, P_FDR := p.adjust(p, method = "BH")]
phenotype_labels <- c(SmokingInitiation = "Smoking initiation", CigarettesPerDay = "Cigarettes per day")
p2s_valid[, Phenotype := phenotype_labels[phenotype]]
p2s_valid[, `Outcome type` := ifelse(outcome_type == "binary", "Binary", "Continuous")]
p2s_valid[, Method := ifelse(n_snps == 1, "Wald ratio", "IVW")]
p2s_valid[, `Effect estimate (95% CI)` := ifelse(outcome_type == "binary",
  sprintf("OR=%.2f (%.2f-%.2f)", OR_or_NA, OR_lo_or_NA, OR_hi_or_NA),
  sprintf("beta=%.4f SD (%.4f-%.4f)", beta, beta_lo_if_continuous, beta_hi_if_continuous))]
p2s_valid[, `WME (95% CI)` := ifelse(is.na(wme_beta), "-", sprintf("%.4f (%.4f-%.4f)", wme_beta, wme_beta-1.96*wme_se, wme_beta+1.96*wme_se))]
p2s_valid[, `MBE (95% CI)` := ifelse(is.na(mbe_beta), "-", sprintf("%.4f (%.4f-%.4f)", mbe_beta, mbe_beta-1.96*mbe_se, mbe_beta+1.96*mbe_se))]

table2 <- p2s_valid[, .(
  Protein = gene, Phenotype, `Outcome type`, Method, nSNP = n_snps, `Effect estimate (95% CI)`,
  P = fmt_p(p), `P (FDR-adjusted)` = fmt_p(P_FDR),
  `Weighted median (95% CI)` = `WME (95% CI)`, `WME P` = fmt_p(wme_p),
  `Weighted mode (95% CI)` = `MBE (95% CI)`, `MBE P` = fmt_p(mbe_p),
  `Heterogeneity P` = ifelse(is.na(heterogeneity_p), "-", sprintf("%.3f", heterogeneity_p)),
  `Egger intercept P` = ifelse(is.na(egger_intercept_p), "-", sprintf("%.3f", egger_intercept_p))
)]
p2s_pheno_order <- c("Smoking initiation", "Cigarettes per day")
table2[, Phenotype := factor(Phenotype, levels = p2s_pheno_order)]
setorder(table2, Protein, Phenotype)
table2[, Phenotype := as.character(Phenotype)]

con <- file(file.path(OUT_DIR, "Supplementary_Table_2.tsv"), "w")
writeLines(c("Supplementary Table 2. MR estimates for 14 panel proteins and smoking behaviour.",
             "Cis-pQTL instruments (as Table 1) vs smoking initiation (n=249,171: 140,856 ever-smokers, 108,315 never-smokers) and cigarettes per day (n=143,210), GSCAN excl. UKB/23andMe. FDR (BH) across all tests.", ""), con)
close(con)
fwrite(table2, file.path(OUT_DIR, "Supplementary_Table_2.tsv"), append = TRUE, col.names = TRUE, sep = "\t")

# ============================================================
# SECTION 6: scope and clump smoking-instrument loci (reverse direction)
# ============================================================

scope_smoking_loci <- function(pheno_file) {
  f <- file.path(DATA_DIR, "smoking_gwas", pheno_file)
  cmd <- paste(DECOMPRESS, shQuote(f), "| awk -F'\t' 'NR>1 && $8<5e-8 {print}'")
  out_lines <- system(cmd, intern = TRUE)
  if (length(out_lines) == 0) return(NULL)
  header <- strsplit(system(paste(DECOMPRESS, shQuote(f), "| head -1"), intern = TRUE), "\t")[[1]]
  d <- fread(text = paste(out_lines, collapse = "\n"), header = FALSE, sep = "\t")
  setnames(d, seq_along(header), header)
  setorder(d, CHROM, POS)
  d
}
smkinit_gwsig <- scope_smoking_loci(smoking_files["SmokingInitiation"])
fwrite(smkinit_gwsig, file.path(OUT_DIR, "smoking_initiation_gwsig_raw.csv"))
cigday_gwsig <- scope_smoking_loci(smoking_files["CigarettesPerDay"])
fwrite(cigday_gwsig, file.path(OUT_DIR, "cigarettes_per_day_gwsig_raw.csv"))

# Locus boundaries determined empirically from the scoping step above
# (gap > 1Mb between consecutive SNPs on a chromosome = separate locus).
smoking_loci <- list(
  SmokingInitiation = list(
    ld_dir = file.path(DATA_DIR, "ld", "smoking_initiation"),
    cand = smkinit_gwsig,
    out_combined = file.path(OUT_DIR, "smoking_initiation_final_instruments.csv"),
    loci = data.table(
      name  = c("locus_chr1","locus_chr2a","locus_chr2b","locus_chr3","locus_chr4a",
                "locus_chr4b","locus_chr5a","locus_chr5b","locus_chr11","locus_chr15"),
      chrom = c("1","2","2","3","4","4","5","5","11","15"),
      start = c(44007648,146111968,155664656,49911155,68020214,140884723,154784177,166987674,111992273,83602849),
      end   = c(44086831,146156679,155723912,50224225,68038323,140939110,154882268,167006676,112973497,83977166)
    )
  ),
  CigarettesPerDay = list(
    ld_dir = file.path(DATA_DIR, "ld", "cigarettes_per_day"),
    cand = cigday_gwsig,
    out_combined = file.path(OUT_DIR, "cigarettes_per_day_final_instruments.csv"),
    loci = data.table(
      name  = c("locus_chr7","locus_chr8","locus_chr9","locus_chr11",
                "locus_chr15","locus_chr16","locus_chr19","locus_chr20"),
      chrom = c("7","8","9","11","15","16","19","20"),
      start = c(32261458,42510343,136502321,16231392,78243579,52074123,41207206,61986949),
      end   = c(32379218,42661847,136505241,16231392,79215568,52125009,41439939,61992005)
    )
  )
)

clump_smoking_locus <- function(locus_name, chrom, start, end, cand_all, ld_dir) {
  cand <- cand_all[as.character(CHROM) == chrom & POS >= start - 10000 & POS <= end + 10000]
  ld_file <- file.path(ld_dir, paste0(locus_name, "_ld.csv"))
  if (!file.exists(ld_file)) return(NULL)
  ld <- read_ld_matrix(ld_file)
  independent_rsids <- greedy_clump(cand$RSID, cand$PVALUE, ld$mat, ld$labels, sig_is_log10 = FALSE)
  final <- cand[RSID %in% independent_rsids]
  fwrite(final, file.path(ld_dir, paste0(locus_name, "_final_instruments.csv")))
  final
}

for (pheno_cfg in smoking_loci) {
  loci <- pheno_cfg$loci
  mapply(function(name, chrom, start, end) clump_smoking_locus(name, chrom, start, end, pheno_cfg$cand, pheno_cfg$ld_dir),
         loci$name, loci$chrom, loci$start, loci$end, SIMPLIFY = FALSE)
  all_final <- rbindlist(lapply(loci$name, function(nm) {
    f <- file.path(pheno_cfg$ld_dir, paste0(nm, "_final_instruments.csv"))
    if (file.exists(f)) { d <- fread(f); d$locus <- nm; d } else NULL
  }), fill = TRUE)
  fwrite(all_final, pheno_cfg$out_combined)
}

# ============================================================
# SECTION 7: extract smoking-instrument effects on protein levels
# ============================================================

smoking_init_exp <- fread(file.path(OUT_DIR, "smoking_initiation_final_instruments.csv")); smoking_init_exp$exposure_phenotype <- "SmokingInitiation"
cigday_exp <- fread(file.path(OUT_DIR, "cigarettes_per_day_final_instruments.csv")); cigday_exp$exposure_phenotype <- "CigarettesPerDay"
smoking_instruments_all <- rbindlist(list(smoking_init_exp, cigday_exp), fill = TRUE)

needed_chroms <- unique(as.character(smoking_instruments_all$CHROM))
snp_lookup <- rbindlist(lapply(needed_chroms, function(chrom) {
  map <- get_snp_map(chrom)
  if (is.null(map)) return(NULL)
  map[rsid %in% smoking_instruments_all[CHROM == chrom, RSID]]
}), fill = TRUE)
target_ids <- snp_lookup$ID

extract_smoking_in_protein <- function(gene) {
  if (!gene %in% names(pqtl_tar_files)) return(NULL)
  contents <- untar(pqtl_tar_files[gene], list = TRUE)
  results <- list()
  for (chrom in needed_chroms) {
    target <- grep(paste0("discovery_chr", chrom, "_"), contents, value = TRUE)
    if (length(target) == 0) next
    tmp_dir <- tempfile(); dir.create(tmp_dir)
    untar(pqtl_tar_files[gene], files = target, exdir = tmp_dir)
    pqtl <- fread(file.path(tmp_dir, target))
    hit <- pqtl[ID %in% target_ids]
    if (nrow(hit) > 0) { hit[, PVAL := 10^(-LOG10P)]; hit$gene <- gene; results[[chrom]] <- hit }
    unlink(tmp_dir, recursive = TRUE)
  }
  if (length(results) == 0) return(NULL)
  rbindlist(results, fill = TRUE)
}

smoking_to_protein_outcome <- rbindlist(lapply(GENES, extract_smoking_in_protein), fill = TRUE)
smoking_to_protein_outcome <- merge(smoking_to_protein_outcome, snp_lookup[, .(ID, rsid)], by = "ID", all.x = TRUE)
smoking_to_protein_outcome <- merge(smoking_to_protein_outcome, smoking_instruments_all[, .(RSID, exposure_phenotype)],
                                     by.x = "rsid", by.y = "RSID", all.x = TRUE)
fwrite(smoking_to_protein_outcome, file.path(OUT_DIR, "smoking_to_protein_outcome_data.csv"))

# ============================================================
# SECTION 8: MR, smoking -> protein (with Radial MR outlier detection)
# ============================================================

mr_smoking_to_protein <- function(target_gene, target_phenotype) {
  exp_raw <- smoking_instruments_all[exposure_phenotype == target_phenotype]
  out_raw <- smoking_to_protein_outcome[smoking_to_protein_outcome$gene == target_gene & smoking_to_protein_outcome$exposure_phenotype == target_phenotype]
  if (nrow(exp_raw) == 0 || nrow(out_raw) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "no data"))

  exp_df <- as.data.frame(exp_raw)
  if (target_phenotype == "SmokingInitiation") {
    exp_df$ncase <- SMKINIT_NCASE; exp_df$ncontrol <- SMKINIT_NCONTROL; exp_df$samplesize <- SMKINIT_NCASE + SMKINIT_NCONTROL
    exp_dat <- TwoSampleMR::format_data(exp_df, type = "exposure",
      snp_col = "RSID", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF",
      eaf_col = "AF", pval_col = "PVALUE", ncase_col = "ncase", ncontrol_col = "ncontrol", samplesize_col = "samplesize")
  } else {
    exp_df$samplesize <- exp_df$N
    exp_dat <- TwoSampleMR::format_data(exp_df, type = "exposure",
      snp_col = "RSID", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALT", other_allele_col = "REF",
      eaf_col = "AF", pval_col = "PVALUE", samplesize_col = "samplesize")
  }
  exp_dat$exposure <- target_phenotype

  out_df <- as.data.frame(out_raw); out_df$samplesize <- UKBPPP_N
  out_dat <- TwoSampleMR::format_data(out_df, type = "outcome",
    snp_col = "rsid", beta_col = "BETA", se_col = "SE", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0",
    eaf_col = "A1FREQ", pval_col = "PVAL", samplesize_col = "samplesize")
  out_dat$outcome <- target_gene

  harm <- tryCatch(TwoSampleMR::harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
  if (is.null(harm) || nrow(harm) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "harmonisation failed"))
  harm <- harm[harm$mr_keep, ]
  if (nrow(harm) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "all SNPs dropped"))

  harm_steiger <- tryCatch(TwoSampleMR::steiger_filtering(harm), error = function(e) NULL)
  if (is.null(harm_steiger)) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "steiger failed"))
  harm_final <- harm_steiger[!is.na(harm_steiger$steiger_dir) & harm_steiger$steiger_dir, ]
  if (nrow(harm_final) == 0) return(data.table(gene = target_gene, phenotype = target_phenotype, status = "failed Steiger"))

  n_snps <- nrow(harm_final)
  mrin <- MendelianRandomization::mr_input(bx = harm_final$beta.exposure, bxse = harm_final$se.exposure,
                                            by = harm_final$beta.outcome, byse = harm_final$se.outcome)
  ivw <- MendelianRandomization::mr_ivw(mrin)
  egger <- NULL; wme <- NULL; mbe <- NULL
  if (n_snps >= 3) {
    egger <- tryCatch(MendelianRandomization::mr_egger(mrin), error = function(e) NULL)
    wme   <- tryCatch(MendelianRandomization::mr_median(mrin, weighting = "weighted"), error = function(e) NULL)
    mbe   <- tryCatch(MendelianRandomization::mr_mbe(mrin), error = function(e) NULL)
  }
  het_p <- if (n_snps == 1) NA else tryCatch(ivw@Heter.Stat[2], error = function(e) NA)

  n_outliers <- 0; outlier_snps <- NA_character_; beta_excl <- NA; se_excl <- NA; p_excl <- NA
  if (n_snps >= 3) {
    radial <- radial_outlier_test(harm_final$SNP, harm_final$beta.exposure, harm_final$beta.outcome, harm_final$se.outcome)
    n_outliers <- sum(radial$is_outlier)
    if (n_outliers > 0) {
      outlier_snps <- paste(radial[is_outlier == TRUE, rsid], collapse = ";")
      keep <- !(harm_final$SNP %in% radial[is_outlier == TRUE, rsid])
      if (sum(keep) >= 1) {
        clean <- run_ivw_manual(harm_final$beta.exposure[keep], harm_final$beta.outcome[keep], harm_final$se.outcome[keep])
        beta_excl <- clean$beta; se_excl <- clean$se; p_excl <- clean$p
      }
    }
  }

  data.table(
    gene = target_gene, phenotype = target_phenotype, status = "ok", n_snps = n_snps,
    beta = get_slot_safe(ivw, "Estimate"), se = get_slot_safe(ivw, "StdError"), p = get_slot_safe(ivw, "Pvalue"),
    heterogeneity_p = het_p, egger_intercept_p = get_slot_safe(egger, "Pvalue.Int"),
    wme_beta = get_slot_safe(wme, "Estimate"), wme_se = get_slot_safe(wme, "StdError"), wme_p = get_slot_safe(wme, "Pvalue"),
    mbe_beta = get_slot_safe(mbe, "Estimate"), mbe_se = get_slot_safe(mbe, "StdError"), mbe_p = get_slot_safe(mbe, "Pvalue"),
    n_radial_outliers = n_outliers, radial_outlier_snps = outlier_snps,
    beta_excl_outliers = beta_excl, se_excl_outliers = se_excl, p_excl_outliers = p_excl
  )
}

s2p_phenotypes <- c("SmokingInitiation", "CigarettesPerDay")
s2p_results <- rbindlist(lapply(GENES, function(g) rbindlist(lapply(s2p_phenotypes, function(p) mr_smoking_to_protein(g, p)), fill = TRUE)), fill = TRUE)
fwrite(s2p_results, file.path(OUT_DIR, "mr_smoking_to_protein.csv"))

# ---- Table 3 (FDR separately per phenotype - distinct causal questions) ----
s2p_valid <- s2p_results[status == "ok"]
s2p_valid[, P_FDR := p.adjust(p, method = "BH"), by = phenotype]
s2p_valid[, Phenotype := phenotype_labels[phenotype]]
s2p_valid[, Method := ifelse(n_snps == 1, "Wald ratio", "IVW")]
s2p_valid[, `Beta (95% CI)` := sprintf("%.4f (%.4f-%.4f)", beta, beta-1.96*se, beta+1.96*se)]
s2p_valid[, `WME (95% CI)` := ifelse(is.na(wme_beta), "-", sprintf("%.4f (%.4f-%.4f)", wme_beta, wme_beta-1.96*wme_se, wme_beta+1.96*wme_se))]
s2p_valid[, `MBE (95% CI)` := ifelse(is.na(mbe_beta), "-", sprintf("%.4f (%.4f-%.4f)", mbe_beta, mbe_beta-1.96*mbe_se, mbe_beta+1.96*mbe_se))]
s2p_valid[, `Radial MR outlier(s)` := ifelse(is.na(radial_outlier_snps) | radial_outlier_snps == "", "None", radial_outlier_snps)]
s2p_valid[, `Beta excl. outlier(s) (95% CI)` := ifelse(is.na(beta_excl_outliers), "-",
  sprintf("%.4f (%.4f-%.4f)", beta_excl_outliers, beta_excl_outliers-1.96*se_excl_outliers, beta_excl_outliers+1.96*se_excl_outliers))]

table3 <- s2p_valid[, .(
  Protein = gene, Phenotype, Method, nSNP = n_snps, `Beta (95% CI)`,
  P = fmt_p(p), `P (FDR-adjusted)` = fmt_p(P_FDR),
  `Weighted median (95% CI)` = `WME (95% CI)`, `WME P` = fmt_p(wme_p),
  `Weighted mode (95% CI)` = `MBE (95% CI)`, `MBE P` = fmt_p(mbe_p),
  `Heterogeneity P` = ifelse(is.na(heterogeneity_p), "-", sprintf("%.3f", heterogeneity_p)),
  `Egger intercept P` = ifelse(is.na(egger_intercept_p), "-", sprintf("%.3f", egger_intercept_p)),
  `Radial MR outlier(s)`, `Beta excl. outlier(s) (95% CI)`, `P excl. outlier(s)` = fmt_p(p_excl_outliers)
)]
table3[, Phenotype := factor(Phenotype, levels = p2s_pheno_order)]
setorder(table3, Protein, Phenotype)
table3[, Phenotype := as.character(Phenotype)]

con <- file(file.path(OUT_DIR, "Supplementary_Table_3.tsv"), "w")
writeLines(c("Supplementary Table 3. MR estimates for the causal effect of smoking initiation and cigarettes per day on plasma levels of 14 candidate proteins.",
             "Genetic instruments for smoking initiation (10 loci) and cigarettes per day (8 loci), GSCAN excl. UKB/23andMe, vs UKB-PPP protein levels (n=34,557). Beta = SD change in protein level; not exponentiated. Radial MR (Bowden et al. 2018) flags outlier instruments via Cochran's Q contribution (Bonferroni-corrected); outlier-excluded estimates shown where detected. FDR (BH) separately per phenotype.", ""), con)
close(con)
fwrite(table3, file.path(OUT_DIR, "Supplementary_Table_3.tsv"), append = TRUE, col.names = TRUE, sep = "\t")

# ============================================================
# SECTION 9: Figure 3 - three-panel forest plot
# ============================================================

GENE_ORDER <- c("CEACAM5","WFDC2","LAMP3","GDF15","CXCL17","ALPP","MMP12",
                "PIGR","PLAUR","PRSS8","SFTPD","CDCP1","SFTPA1","TNFSF13B")

fmt_p_hybrid <- function(p) {
  s <- formatC(p, format = "e", digits = 2)
  parts <- strsplit(s, "e")
  mantissa <- sapply(parts, `[`, 1); exponent <- as.integer(sapply(parts, `[`, 2))
  ifelse(p < 0.01, sprintf("%s %%*%% 10^%d", mantissa, exponent), sprintf("%.3f", p))
}

panel_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
        axis.ticks = element_line(colour = "black"), plot.title = element_text(size = 12, face = "bold"),
        aspect.ratio = 1, legend.position = "none")

d_a <- lc_results_df[status == "ok"]
d_a[, P_FDR := p.adjust(p, method = "BH")]
d_a <- d_a[outcome == "lung_cancer"]
d_a[, sig := P_FDR < 0.05]
d_a[, gene := factor(gene, levels = rev(GENE_ORDER))]
d_a[, label_expr := sprintf("italic(P)[FDR] == %s", fmt_p_hybrid(P_FDR))]
label_x_a <- max(d_a$OR_hi, na.rm = TRUE) * 1.05

panel_a <- ggplot(d_a, aes(x = OR, y = gene)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "black") +
  geom_errorbarh(aes(xmin = OR_lo, xmax = OR_hi, colour = sig), height = 0.25, linewidth = 0.55) +
  geom_point(aes(colour = sig), size = 2.8) +
  geom_text(aes(x = label_x_a, label = label_expr), parse = TRUE, hjust = 0, size = 2.9, colour = "black") +
  scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "black"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.6))) +
  labs(x = "Odds ratio (per 1-SD increase in protein level)", y = NULL,
       title = "Causal effect of the Pandya et al. 14-protein signature on lung cancer risk") +
  panel_theme + theme(plot.margin = margin(5.5, 20, 5.5, 5.5))

d_bc <- s2p_results[status == "ok"]
d_bc[, P_FDR := p.adjust(p, method = "BH"), by = phenotype]
d_bc[, sig := P_FDR < 0.05]
d_bc[, ci_lo := beta - 1.96*se]; d_bc[, ci_hi := beta + 1.96*se]
d_bc[, gene := factor(gene, levels = rev(GENE_ORDER))]
d_bc[, label_expr := sprintf("italic(P)[FDR] == %s", fmt_p_hybrid(P_FDR))]

make_forest_bc <- function(pheno_key, panel_title, x_lab) {
  sub <- d_bc[phenotype == pheno_key]
  label_x <- max(sub$ci_hi) * 1.05
  ggplot(sub, aes(x = beta, y = gene)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "black") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, colour = sig), height = 0.25, linewidth = 0.55) +
    geom_point(aes(colour = sig), size = 2.8) +
    geom_text(aes(x = label_x, label = label_expr), parse = TRUE, hjust = 0, size = 2.9, colour = "black") +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "black"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.6))) +
    labs(x = x_lab, y = NULL, title = panel_title) +
    panel_theme + theme(plot.margin = margin(5.5, 20, 5.5, 20))
}

panel_b <- make_forest_bc("SmokingInitiation", "Genetic effect of smoking initiation on the Pandya et al. 14-protein signature",
                           "SNP effect on protein level (SD units) per log-odds of smoking initiation")
panel_c <- make_forest_bc("CigarettesPerDay", "Genetic effect of cigarettes per day on the Pandya et al. 14-protein signature",
                           "SNP effect on protein level (SD units) per SD of cigarettes per day")

combined <- (panel_a | panel_b | panel_c) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT_DIR, "figure3_three_panel.pdf"), combined, width = 20, height = 7, units = "in")
ggsave(file.path(OUT_DIR, "figure3_three_panel.png"), combined, width = 20, height = 7, units = "in", dpi = 200)
