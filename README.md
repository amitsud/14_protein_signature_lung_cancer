# Reanalysis code: correspondence on Pandya et al. (Cell 2026)

Reproducible R code accompanying a correspondence on Pandya et al., "Plasma signals of lung tumor promotion for molecular cancer prevention," Cell 2026; DOI 10.1016/j.cell.2026.05.005.

The scripts reproduce the statistical reanalyses and figures in the correspondence: predictive/screening conversions, the CANTOS treatment-interaction test, and the disease-association figures derived from the UK Biobank plasma proteome atlas.

## Scripts

`detection_rate_conversions.R` converts reported ROC-AUCs into screening terms — detection rate at a 5% false-positive rate (DR5) and post-test odds — under an equal-variance binormal ROC model.

`interaction_test_altman_bland.R` tests whether the canakinumab effect differs between high- and low-signature CANTOS subgroups, using the Altman–Bland method for the difference between two estimates (ratio of odds ratios, with z and p).

`protein_forest_plot.R` draws a per-protein forest plot of the 14 signature proteins' hazard ratios across six endpoints (lung cancer, smoking dependency, early- and late-onset COPD, idiopathic pulmonary fibrosis, peripheral vascular disease).

`barcode_enrichment_plot.R` draws a rank-enrichment ("barcode") plot showing where the 14 proteins rank among all plasma proteins as risk markers for each disease, with a Wilcoxon/Mann–Whitney AUC (Hanley–McNeil 95% CI) and a permutation test.

The first two scripts are self-contained (they take published summary statistics as inputs and need no external data). The two figure scripts require the atlas data files described below.

## Data

The figure scripts read per-protein, smoking-adjusted incident-disease hazard ratios from the UK Biobank Proteome–Phenome Atlas (Deng et al., Cell 2024; DOI 10.1016/j.cell.2024.10.045), available via the portal at https://proteome-phenome-atlas.com/. The underlying UK Biobank resource is available to approved researchers at https://www.ukbiobank.ac.uk/.

These data are not redistributed here. To run the figure scripts, export the relevant per-disease tables from the atlas and place them in `data/ukbiobank/`, one file per disease:

    data/ukbiobank/
      lungcancer.txt  earlycopd.txt  latecopd.txt  ild.txt  smokingdependency.txt
      ihd.txt  heartfailure.txt  ischaemicva.txt  pvd.txt

Each file is tab-delimited with columns including `Protein`, `Disease`, `NB_case`, `HR[95%CI]` (e.g. `1.58 [1.32-1.91]`) and `P_value`. The scripts parse these automatically.

## Running

With R installed, from the repository root:

    install.packages(c("readr","dplyr","tidyr","stringr","ggplot2","purrr","tibble","forcats"))

    Rscript detection_rate_conversions.R
    Rscript interaction_test_altman_bland.R
    Rscript protein_forest_plot.R          # needs data/ukbiobank/
    Rscript barcode_enrichment_plot.R      # needs data/ukbiobank/

Figures and tables are written to `outputs/`. The figure scripts print the disease endpoint they selected for each file — check this matches the intended endpoint.

## Notes and caveats

DR5 values are derived from reported AUCs under an equal-variance binormal ROC model; they are not values reported by the original authors. The Altman–Bland test uses standard errors back-calculated from the published subgroup confidence intervals. In the barcode analysis the 14 proteins are correlated, so the permutation p-values and AUC confidence intervals are approximate (anti-conservative). Lung cancer, COPD, idiopathic pulmonary fibrosis and smoking endpoints are not independent of the signature's derivation and are shown for reference.

## Citation

If you use this code, please cite the correspondence (details to follow on publication) and the underlying data sources above.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The data are not included and remain subject to the terms of the UK Biobank Proteome–Phenome Atlas and UK Biobank.
