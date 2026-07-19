# Reanalysis code: correspondence on Pandya et al. (Cell 2026)

Reproducible R code accompanying a correspondence on Pandya et al., "Plasma signals of lung tumor promotion for molecular cancer prevention," Cell 2026; DOI 10.1016/j.cell.2026.05.005.

The scripts reproduce the statistical reanalyses and figures in the correspondence: predictive/screening conversions, the CANTOS treatment-interaction test, the disease-association figures derived from the UK Biobank plasma proteome atlas, and the Mendelian randomisation (MR) analyses testing the proposed causal mechanism.

## Scripts

`detection_rate_conversions.R` converts reported ROC-AUCs into screening terms — detection rate at a 5% false-positive rate (DR5) and post-test odds — under an equal-variance binormal ROC model.

`interaction_test_altman_bland.R` tests whether the canakinumab effect differs between high- and low-signature CANTOS subgroups, using the Altman–Bland method for the difference between two estimates (ratio of odds ratios, with z and p).

`protein_forest_plot.R` draws a per-protein forest plot of the 14 signature proteins' hazard ratios across six endpoints (lung cancer, smoking dependency, early- and late-onset COPD, idiopathic pulmonary fibrosis, peripheral vascular disease).

`barcode_enrichment_plot.R` draws a rank-enrichment ("barcode") plot showing where the 14 proteins rank among all plasma proteins as risk markers for each disease, with a Wilcoxon/Mann–Whitney AUC (Hanley–McNeil 95% CI) and a permutation test.

`mendelian_randomization_pipeline.R` runs the two-sample MR analyses testing the correspondence's causal claim: whether the 14 panel proteins causally affect lung cancer risk or smoking behaviour, and whether smoking behaviour causally affects the panel proteins. Uses cis-pQTLs from UKB-PPP as protein instruments and independent, LD-clumped genome-wide-significant loci from GSCAN as smoking instruments (see "LD clumping" below). Produces Supplementary Tables 1–3 (protein → lung cancer; protein → smoking; smoking → protein) and Figure 3 (three-panel forest plot).

The first two scripts are self-contained (they take published summary statistics as inputs and need no external data). The two figure scripts and the MR pipeline require the data files described below.

## Data

### UK Biobank Proteome–Phenome Atlas

`protein_forest_plot.R` and `barcode_enrichment_plot.R` read per-protein, smoking-adjusted incident-disease hazard ratios from the UK Biobank Proteome–Phenome Atlas (Deng et al., Cell 2024; DOI 10.1016/j.cell.2024.10.045), available via the portal at https://proteome-phenome-atlas.com/. The underlying UK Biobank resource is available to approved researchers at https://www.ukbiobank.ac.uk/.

These data are not redistributed here. To run the figure scripts, export the relevant per-disease tables from the atlas and place them in `data/ukbiobank/`, one file per disease:

    data/ukbiobank/
      lungcancer.txt  earlycopd.txt  latecopd.txt  ild.txt  smokingdependency.txt
      ihd.txt  heartfailure.txt  ischaemicva.txt  pvd.txt

Each file is tab-delimited with columns including `Protein`, `Disease`, `NB_case`, `HR[95%CI]` (e.g. `1.58 [1.32-1.91]`) and `P_value`. The scripts parse these automatically.

### MR pipeline inputs

`mendelian_randomization_pipeline.R` requires several raw GWAS resources, none of which are redistributed here, placed under `data/`:

    data/
      pqtl_raw/                          UKB-PPP (Sun et al. 2023) per-protein .tar files, one per panel protein
      snp_map/                           UKB-PPP olink_rsid_map_*.tsv.gz files (one per chromosome)
      smoking_gwas/
        SmokingInitiation.WithoutUKB.txt.gz
        CigarettesPerDay.WithoutUKB.txt.gz   (GSCAN, UK Biobank/23andMe excluded)
      lung_gwas/
        lung_cancer/harmonized/            ILCCO/TRICL summary statistics (GWAS-SSF harmonised format)
        lung_adenocarcinoma/harmonized/
        lung_squamous/harmonized/
        lung_eversmoker/harmonized/
        lung_neversmoker/harmonized/
      ld/
        protein/                           per-protein LD correlation matrices (PolyFun/Weissbrod UK Biobank EUR panel)
        smoking_initiation/                 per-locus LD matrices for the 10 Smoking Initiation loci
        cigarettes_per_day/                 per-locus LD matrices for the 8 Cigarettes Per Day loci

UKB-PPP data require UK Biobank access (as above). GSCAN summary statistics (excluding UK Biobank and 23andMe) are available from the GWAS & Sequencing Consortium of Alcohol and Nicotine Use. ILCCO/TRICL lung cancer summary statistics are available via the GWAS Catalog. LD matrices were extracted from the PolyFun/Weissbrod UK Biobank European-ancestry reference panel (337,545 individuals); this project's own extraction step is not included here (see caveats).

## Running

With R installed, from the repository root:

    install.packages(c("readr","dplyr","tidyr","stringr","ggplot2","purrr","tibble","forcats",
                        "data.table","remotes","MendelianRandomization","patchwork"))
    remotes::install_github("MRCIEU/TwoSampleMR")

    Rscript detection_rate_conversions.R
    Rscript interaction_test_altman_bland.R
    Rscript protein_forest_plot.R              # needs data/ukbiobank/
    Rscript barcode_enrichment_plot.R           # needs data/ukbiobank/
    Rscript mendelian_randomization_pipeline.R  # needs data/ as described above; can take a while (large GWAS files)

Figures and tables from the first four scripts are written to `outputs/`. The MR pipeline writes its outputs to `results/`. The figure scripts print the disease endpoint they selected for each file — check this matches the intended endpoint.

## Notes and caveats

DR5 values are derived from reported AUCs under an equal-variance binormal ROC model; they are not values reported by the original authors. The Altman–Bland test uses standard errors back-calculated from the published subgroup confidence intervals. In the barcode analysis the 14 proteins are correlated, so the permutation p-values and AUC confidence intervals are approximate (anti-conservative). Lung cancer, COPD, idiopathic pulmonary fibrosis and smoking endpoints are not independent of the signature's derivation and are shown for reference.

The MR pipeline's LD-matrix extraction step is not included in this repository; the pipeline expects pre-extracted per-locus LD matrices as described above. 

## Citation

If you use this code, please cite the correspondence (details to follow on publication) and the underlying data sources above.

## License

Code in this repository is released under the MIT License (see `LICENSE`). The data are not included and remain subject to the terms of the UK Biobank Proteome–Phenome Atlas, UK Biobank, GSCAN, and the GWAS Catalog.
