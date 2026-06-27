################################################################################
## Barcode enrichment of the 14-protein signature.
## For each disease, all proteins are ranked by association strength (z=log(HR)/SE);
## each signature protein is a tick at its rank percentile (1 = strongest, left).
## Ticks bunched left = the signature is enriched among that disease's top risk
## markers. AUC = probability a signature protein outranks a random protein
## (Wilcoxon / Mann-Whitney), with a Hanley-McNeil 95% CI.
##
## NOTE: figure shows breadth across lung, smoking and cardiovascular disease.
## AUC CI is approximate (proteins are correlated). HRs are smoking-adjusted.
##
## USAGE: place the per-protein atlas .txt files in DATA_DIR (see README), then
## run from the repository root:  Rscript fig_barcode_enrichment.R
## Outputs (figure + table) are written to OUT_DIR.
################################################################################

## ----------------------------- CONFIG -------------------------------------- ##
DATA_DIR <- "data/ukbiobank"   # folder containing the per-protein atlas .txt files
OUT_DIR  <- "outputs"          # folder for figures and tables
B_PERM   <- 5000
ADJ <- "adjusted for age, sex, ethnicity, Townsend index, BMI, smoking status, fasting time, season and blood age"

SIGNATURE <- toupper(c("SFTPA1","SFTPD","CXCL17","WFDC2","LAMP3","PIGR","PRSS8",
                       "GDF15","PLAUR","CDCP1","TNFSF13B","CEACAM5","ALPP","MMP12"))
alias <- c("CXL17"="CXCL17","U-PAR"="PLAUR","UPAR"="PLAUR")

diseases <- tibble::tribble(
  ~file,                  ~label,                          ~keyword,                        ~group,
  "lungcancer.txt",       "Non-small cell lung cancer",    "non-small cell|nsclc",          "Lung",
  "earlycopd.txt",        "COPD (early onset)",            "copd",                          "Lung",
  "latecopd.txt",         "COPD (late onset)",             "copd",                          "Lung",
  "ild.txt",              "Idiopathic pulmonary fibrosis", "idiopathic|fibros|ipf",         "Lung",
  "smokingdependency.txt","Smoking dependency",            "smok|tobacco|nicotine|depend",  "Smoking",
  "ihd.txt",              "Ischaemic heart disease",       "ischaem.*heart|coronary|heart", "Cardiovascular",
  "heartfailure.txt",     "Heart failure",                 "heart failure|cardiac failure", "Cardiovascular",
  "ischaemicva.txt",      "Ischaemic stroke",              "stroke|cerebrovasc",            "Cardiovascular",
  "pvd.txt",              "Peripheral vascular disease",   "periph|vascular",               "Cardiovascular"
)

## --------------------------- LIBRARIES ------------------------------------- ##
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(purrr); library(forcats); library(tibble)
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
set.seed(1)

theme_letter <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour = "grey92"),
        legend.position = "top")

## --------------------------- HELPERS --------------------------------------- ##
parse_hr <- function(s) {
  s <- as.character(s); s <- gsub("\u2013|\u2212", "-", s)
  m <- str_match(s, "([0-9]*\\.?[0-9]+)\\s*\\[\\s*([0-9]*\\.?[0-9]+)\\s*-\\s*([0-9]*\\.?[0-9]+)\\s*\\]")
  tibble(hr = as.numeric(m[,2]), ci_low = as.numeric(m[,3]), ci_high = as.numeric(m[,4]))
}
standardize_cols <- function(df) {
  k <- tolower(names(df))
  for (i in seq_along(k)) {
    if (grepl("definition|protein_def", k[i]))        names(df)[i] <- "protein_def"
    else if (k[i] == "protein")                       names(df)[i] <- "protein"
    else if (grepl("category|disease_cat", k[i]))     names(df)[i] <- "disease_category"
    else if (k[i] == "disease")                       names(df)[i] <- "disease"
    else if (grepl("individual|nb_indiv", k[i]))      names(df)[i] <- "n_individual"
    else if (grepl("case", k[i]))                     names(df)[i] <- "n_case"
    else if (grepl("hr", k[i]) && grepl("ci", k[i]))  names(df)[i] <- "hr_string"
    else if (grepl("^hr", k[i]))                      names(df)[i] <- "hr_string"
    else if (grepl("p.?val|^p$|p_value", k[i]))       names(df)[i] <- "p_value"
  }
  df
}
read_atlas <- function(path) {
  if (!file.exists(path)) { warning("missing: ", basename(path)); return(NULL) }
  df <- suppressMessages(read_tsv(path, col_types = cols(.default = "c"),
                                  locale = locale(encoding = "latin1"), progress = FALSE))
  df <- standardize_cols(df); df <- bind_cols(df, parse_hr(df$hr_string))
  df %>% mutate(protein = toupper(trimws(protein)),
                protein = ifelse(protein %in% names(alias), alias[protein], protein),
                disease = trimws(disease),
                n_case  = suppressWarnings(as.numeric(gsub(",", "", n_case))),
                se = (log(ci_high) - log(ci_low)) / (2 * 1.96),
                z  = ifelse(is.finite(se) & se > 0, log(hr) / se, NA_real_))
}
pick_endpoint <- function(d, keyword) {
  labs <- unique(d$disease); if (length(labs) == 1) return(labs)
  if (!is.na(keyword)) {
    hit <- labs[grepl(keyword, labs, ignore.case = TRUE)]
    if (length(hit)) { sub <- d %>% filter(disease %in% hit); return(sub$disease[which.max(sub$n_case)]) }
  }
  ag <- d %>% group_by(disease) %>% summarise(mc = max(n_case, na.rm = TRUE), .groups = "drop")
  ag$disease[which.max(ag$mc)]
}
# Hanley-McNeil 95% CI for an AUC (Mann-Whitney). Approximate (assumes independence).
auc_ci <- function(auc, n1, n2) {
  q1 <- auc / (2 - auc); q2 <- 2 * auc^2 / (1 + auc)
  se <- sqrt((auc*(1-auc) + (n1-1)*(q1 - auc^2) + (n2-1)*(q2 - auc^2)) / (n1 * n2))
  c(lo = max(0, auc - 1.96*se), hi = min(1, auc + 1.96*se))
}

## --------------------------- COMPUTE --------------------------------------- ##
cat("\nEndpoints selected (verify):\n")
per_disease <- pmap(diseases, function(file, label, keyword, group) {
  d <- read_atlas(file.path(DATA_DIR, file)); if (is.null(d)) return(NULL)
  used <- pick_endpoint(d, keyword)
  cat(sprintf("  %-30s endpoint: '%s'\n", label, used))
  d <- d %>% filter(disease == used, is.finite(z)) %>% distinct(protein, .keep_all = TRUE)
  N <- nrow(d)
  d$top_pct <- 100 * rank(-d$z, ties.method = "average") / N
  is_sig <- d$protein %in% SIGNATURE; k <- sum(is_sig); if (k < 3) return(NULL)
  ticks <- d %>% filter(is_sig) %>% transmute(label, group, top_pct)
  wt  <- suppressWarnings(wilcox.test(d$z[is_sig], d$z[!is_sig], alternative = "greater"))
  auc <- as.numeric(wt$statistic) / (k * (N - k))
  ci  <- auc_ci(auc, k, N - k)
  obs  <- median(d$top_pct[is_sig])
  null <- replicate(B_PERM, median(sample(d$top_pct, k)))
  list(ticks = ticks,
       stat = tibble(label, group, N_proteins = N, k = k, median_top_pct = obs,
                     AUC = auc, AUC_lo = ci["lo"], AUC_hi = ci["hi"],
                     perm_p = (1 + sum(null <= obs)) / (B_PERM + 1)))
}) %>% compact()

## Clear failure if the data folder is not set up correctly.
if (length(per_disease) == 0)
  stop("No disease files found in '", DATA_DIR,
       "'. Set DATA_DIR (or run from the repo root) so it points to the folder ",
       "containing the atlas .txt files. See README.")

ticks_df <- map_dfr(per_disease, "ticks")
stat_df  <- map_dfr(per_disease, "stat") %>% mutate(perm_p_adj = p.adjust(perm_p, method = "BH"))

cat("\nEnrichment:\n")
print(stat_df %>% arrange(median_top_pct) %>%
        mutate(across(c(AUC, AUC_lo, AUC_hi), ~round(.x,3)),
               median_top_pct = round(median_top_pct,1), perm_p_adj = signif(perm_p_adj,2)) %>%
        select(label, group, AUC, AUC_lo, AUC_hi, median_top_pct, perm_p_adj), n = Inf, width = Inf)
write_csv(stat_df, file.path(OUT_DIR, "table_barcode_enrichment.csv"))

## --------------------------- PLOT ------------------------------------------ ##
ord <- stat_df %>% arrange(desc(median_top_pct)) %>% pull(label)   # most enriched at top
ticks_df <- ticks_df %>% mutate(label = factor(label, levels = ord))
med_df   <- stat_df  %>% mutate(label = factor(label, levels = ord),
                                auclab = sprintf("AUC %.2f (%.2f\u2013%.2f)", AUC, AUC_lo, AUC_hi))

grp_cols <- c("Cardiovascular" = "#d7191c", "Lung" = "#1a9850", "Smoking" = "#2c7bb6")

p <- ggplot(ticks_df, aes(x = top_pct, y = label)) +
  annotate("rect", xmin = 0, xmax = 10, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_vline(xintercept = 50, linetype = 2, colour = "grey70") +
  geom_point(aes(colour = group), shape = 124, size = 5, stroke = 1.1) +
  geom_text(data = med_df, aes(x = 100, y = label, label = auclab), hjust = 1, size = 2.9, colour = "grey20") +
  scale_colour_manual(values = grp_cols, name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0,10,25,50,75,100),
                     labels = c("1\n(top)","10","25","50","75","100\n(bottom)")) +
  labs(x = "Rank percentile among all ~2,900 proteins (1 = strongest risk marker)", y = NULL) +
  theme_letter

ggsave(file.path(OUT_DIR, "fig_barcode_enrichment.pdf"), p, width = 11, height = 5.5, units = "in")
ggsave(file.path(OUT_DIR, "fig_barcode_enrichment.png"), p, width = 11, height = 5.5, units = "in", dpi = 200)

cat("\nFigure: fig_barcode_enrichment.pdf\nTable:  table_barcode_enrichment.csv\n")
