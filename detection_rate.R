################################################################################
## detection_rate_conversions.R
## Convert ROC-AUC into screening terms used in the letter:
##   detection rate at 5% false-positive rate (DR5) and post-test odds.
##
## DR5 from AUC uses the equal-variance binormal ROC model (cases and controls
## Normal with equal variance), so the ROC curve is set by one separation d:
##   AUC = Phi(d/sqrt(2))  =>  d = sqrt(2)*qnorm(AUC);  DR = Phi(d - qnorm(1-fpr))
## DR5 values are derived from reported AUCs under this model, not values the
## original authors reported.
################################################################################

## AUC -> detection rate at false-positive rate `fpr` (default 5%).
auc_to_dr <- function(auc, fpr = 0.05) {
  stopifnot(all(auc > 0.5 & auc < 1))
  d <- sqrt(2) * qnorm(auc)
  pnorm(d - qnorm(1 - fpr))
}

## Post-test odds of disease given a positive result, at prevalence `prev`.
##   OAPR = (DR/FPR) * prev/(1-prev)
post_test <- function(dr, fpr = 0.05, prev = 0.01) {
  post_odds <- (dr / fpr) * (prev / (1 - prev))
  list(post_odds = post_odds,
       ratio_label = sprintf("1:%.1f", 1 / post_odds),
       ppv = post_odds / (1 + post_odds))
}

## ---- Values used in the letter (edit to match cited figures) -------------- ##
auc_inputs <- c(
  "LLPv3 (Fig. 1C)"                  = 0.806,
  "Combined model (Fig. 1C)"        = 0.865,
  "ALPP (single protein, Xiao et al smoking)"  = 0.88,
  "CXCL17 (single protein, Xiao et al smoking)" = 0.87
)

cat("=== AUC -> DR5 (equal-variance binormal approximation) ===\n")
for (nm in names(auc_inputs)) {
  a <- auc_inputs[[nm]]
  cat(sprintf("  %-34s AUC %.3f  ->  DR5 %4.1f%%\n", nm, a, 100 * auc_to_dr(a)))
}

cat("\n=== Post-test odds at 1% prevalence (DR5, 5% FPR) ===\n")
for (nm in c("LLPv3 (Fig. 1C)", "Combined model (Fig. 1C)")) {
  dr <- auc_to_dr(auc_inputs[[nm]])
  pt <- post_test(dr, fpr = 0.05, prev = 0.01)
  cat(sprintf("  %-28s DR5 %4.1f%%  ->  post-test odds %s  (PPV %.1f%%)\n",
              nm, 100 * dr, pt$ratio_label, 100 * pt$ppv))
}
