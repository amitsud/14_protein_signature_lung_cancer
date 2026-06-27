################################################################################
## Per-protein multi-trait forest (6 traits), wide trait spacing, no star labels.
## Each of the 14 signature proteins with its HR for SIX endpoints, showing the
## lung-cancer association is no larger than its associations with smoking,
## airway, fibrotic and systemic vascular disease.
##
## Traits: NSCLC, smoking dependency, COPD (early onset), COPD (late onset),
##         idiopathic pulmonary fibrosis, peripheral vascular disease.
##
## Significance shown by filled (CI excludes 1) vs open (ns) dots. Reported
## P-values are kept in the output table for reference but not drawn.
##
## Spacing knobs:
##   ROWGAP - distance between protein rows (bigger = more space BETWEEN proteins)
##   SPREAD - half-spread of points within a protein (bigger = more space
##            BETWEEN diseases within a protein)
## Source: UK Biobank Proteome-Phenome Atlas; smoking-adjusted incident HRs/SD.
##
## USAGE: place the per-protein atlas .txt files in DATA_DIR (see README), then
## `Rscript fig_per_protein_traits_forest.R`. Outputs are written to OUT_DIR.
################################################################################

## ----------------------------- CONFIG -------------------------------------- ##
DATA_DIR <- "data/ukbiobank"   # folder containing the per-protein atlas .txt files
OUT_DIR  <- "outputs"          # folder for figures and tables
X_BREAKS <- c(0.5, 1, 2, 4, 8)
ROWGAP   <- 3.4    # space between proteins (raise for more)
SPREAD   <- 0.92   # space between diseases within a protein (raise for more)
FIG_H    <- 20     # figure height in inches
ADJ <- "adjusted for age, sex, ethnicity, Townsend index, BMI, smoking status, fasting time, season and blood age"

SIGNATURE <- c("SFTPA1","SFTPD","CXCL17","WFDC2","LAMP3","PIGR","PRSS8",
               "GDF15","PLAUR","CDCP1","TNFSF13B","CEACAM5","ALPP","MMP12")
alias <- c("CXL17"="CXCL17","U-PAR"="PLAUR","UPAR"="PLAUR")

traits <- tibble::tribble(
  ~key,    ~file,                   ~label,                          ~keyword,                       ~color,
  "nsclc", "lungcancer.txt",        "Non-small cell lung cancer",    "non-small cell|nsclc",         "#d7191c",
  "smoke", "smokingdependency.txt", "Smoking dependency",            "smok|tobacco|nicotine|depend", "#2c7bb6",
  "ecopd", "earlycopd.txt",         "COPD (early onset)",            "copd",                         "#1a9850",
  "lcopd", "latecopd.txt",          "COPD (late onset)",             "copd",                         "#66c2a5",
  "ipf",   "ild.txt",               "Idiopathic pulmonary fibrosis", "idiopathic|fibros|ipf",        "#984ea3",
  "pvd",   "pvd.txt",               "Peripheral vascular disease",   "periph|vascular",              "#ff7f00"
)

## --------------------------- LIBRARIES ------------------------------------- ##
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(purrr); library(tibble)
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

theme_letter <- theme_minimal(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
        axis.text.y = element_text(size = 9))

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
  if (!"p_value" %in% names(df)) df$p_value <- NA_character_
  df %>% mutate(protein = toupper(trimws(protein)),
                protein = ifelse(protein %in% names(alias), alias[protein], protein),
                disease = trimws(disease),
                n_case  = suppressWarnings(as.numeric(gsub(",", "", n_case))),
                p       = suppressWarnings(as.numeric(p_value)))
}

## --------------------------- BUILD DATA ------------------------------------ ##
SIG <- toupper(SIGNATURE)
get_trait <- function(file, label, keyword, color) {
  d <- read_atlas(file.path(DATA_DIR, file)); if (is.null(d)) return(NULL)
  d <- d %>% filter(protein %in% SIG, is.finite(hr))
  labs <- unique(d$disease); used <- NA_character_
  if (length(labs) > 1 && !is.na(keyword)) {
    hit <- labs[grepl(keyword, labs, ignore.case = TRUE)]
    if (length(hit)) { sub <- d %>% filter(disease %in% hit); used <- sub$disease[which.max(sub$n_case)] }
  }
  if (is.na(used)) used <- if (length(labs) == 1) labs[1] else {
    ag <- d %>% group_by(disease) %>% summarise(mc = max(n_case, na.rm = TRUE), .groups = "drop")
    ag$disease[which.max(ag$mc)] }
  cat(sprintf("  %-30s endpoint: '%s'%s\n", label, used,
              if (length(labs) > 1) sprintf("  (of %d)", length(labs)) else ""))
  d %>% filter(disease == used) %>%
    transmute(trait = label, protein, hr, ci_low, ci_high, p,
              sig = (ci_low > 1) | (ci_high < 1))
}
cat("\nEndpoints selected (verify):\n")
dat <- pmap_dfr(traits, function(key,file,label,keyword,color) get_trait(file,label,keyword,color)) %>%
  mutate(trait = factor(trait, levels = traits$label))

# order proteins by NSCLC HR (lung-strong at top)
ord <- dat %>% filter(trait == "Non-small cell lung cancer") %>% arrange(hr) %>% pull(protein)
ord <- c(ord, setdiff(SIG, ord))
dat <- dat %>% mutate(protein = factor(protein, levels = ord), idx = as.integer(protein))

# y geometry: protein centres ROWGAP apart; diseases spread +/-SPREAD within
nt <- nrow(traits)
offset  <- setNames(seq(-SPREAD, SPREAD, length.out = nt), traits$label)
dat <- dat %>% mutate(center = idx * ROWGAP, ypos = center + offset[as.character(trait)])

centres <- sort(unique(dat$center))
labs    <- levels(dat$protein)
seps    <- (seq_along(labs)[-length(labs)] + 0.5) * ROWGAP    # dashed lines between proteins
pal     <- setNames(traits$color, traits$label)
cap_h   <- 0.07 * ROWGAP                                       # CI end-cap height (nonzero!)

## ------------------------------ PLOT --------------------------------------- ##
p <- ggplot() +
  geom_hline(yintercept = seps, linetype = 2, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55") +
  geom_errorbarh(data = dat,
                 aes(y = ypos, xmin = ci_low, xmax = ci_high, colour = trait),
                 height = cap_h, linewidth = 0.5, alpha = 0.8) +
  geom_point(data = dat, aes(x = hr, y = ypos, colour = trait, shape = sig),
             size = 2.3, fill = "white", stroke = 0.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21), guide = "none") +
  scale_colour_manual(values = pal, breaks = traits$label) +
  scale_y_continuous(breaks = centres, labels = labs, expand = expansion(add = ROWGAP * 0.6)) +
  scale_x_log10(breaks = X_BREAKS, expand = expansion(mult = c(0.02, 0.06))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(shape = 16, size = 2.8))) +
  labs(x = "Hazard ratio per SD (log scale)", y = NULL) +
  theme_letter

ggsave(file.path(OUT_DIR, "fig_per_protein_traits_forest.pdf"), p, width = 9.5, height = FIG_H, units = "in", limitsize = FALSE)
ggsave(file.path(OUT_DIR, "fig_per_protein_traits_forest.png"), p, width = 9.5, height = FIG_H, units = "in", dpi = 200, limitsize = FALSE)

write_csv(dat %>% select(protein, trait, hr, ci_low, ci_high, p, sig) %>% arrange(protein, trait),
          file.path(OUT_DIR, "table_per_protein_traits.csv"))
cat("\nFigure: fig_per_protein_traits_forest.pdf\nTable:  table_per_protein_traits.csv\n")
