################################################################################
## interaction_test_altman_bland.R
## Test whether two independent effect estimates differ (Altman & Bland, BMJ 2003,
## "Interaction revisited: the difference between two estimates").
##
## For two odds ratios with 95% CIs, back-calculate each log-OR standard error
## from its CI, then test the difference of log-ORs:
##   SE = (log(upper) - log(lower)) / (2 * 1.96)
##   diff = log(OR1) - log(OR2);  SE_diff = sqrt(SE1^2 + SE2^2)
##   z = diff / SE_diff;  p = 2 * pnorm(-|z|)
## exp(diff) is the ratio of odds ratios (ROR), i.e. the interaction effect.
## A non-significant p means the two estimates are not shown to differ; it does
## not show they are equal.
################################################################################

## Compare two ORs (each with 95% CI). Returns ROR, its CI, z and p.
interaction_or <- function(or1, lo1, hi1, or2, lo2, hi2) {
  se1  <- (log(hi1) - log(lo1)) / (2 * 1.96)
  se2  <- (log(hi2) - log(lo2)) / (2 * 1.96)
  diff <- log(or1) - log(or2)
  sed  <- sqrt(se1^2 + se2^2)
  z    <- diff / sed
  list(ROR     = exp(diff),
       ROR_lo  = exp(diff - 1.96 * sed),
       ROR_hi  = exp(diff + 1.96 * sed),
       z       = z,
       p       = 2 * pnorm(-abs(z)))
}

## ---- CANTOS: canakinumab effect, high- vs low-signature subgroups --------- ##
## High-signature OR 0.52 [0.31-0.86]; low-signature OR 0.91 [0.34-2.48].
res <- interaction_or(or1 = 0.52, lo1 = 0.31, hi1 = 0.86,   # high signature
                      or2 = 0.91, lo2 = 0.34, hi2 = 2.48)   # low signature

cat("=== Biomarker-by-treatment interaction (Altman-Bland) ===\n")
cat(sprintf("  Ratio of odds ratios (high vs low): %.2f  [%.2f-%.2f]\n",
            res$ROR, res$ROR_lo, res$ROR_hi))
cat(sprintf("  z = %.2f,  p = %.2f\n", res$z, res$p))
cat("\nThe interaction is not significant: differential benefit of IL-1b blockade\n")
cat("by signature status is not established by these data.\n")
