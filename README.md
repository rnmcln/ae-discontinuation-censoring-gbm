# Loss of survival follow-up after adverse-event-related treatment discontinuation in glioblastoma trials

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24534427.svg)](https://doi.org/10.5281/zenodo.24534427)

Simulation code and review data for a study of the bias that arises when survival follow-up in a single-arm trial ends at treatment discontinuation for adverse events (AEs), while a historical benchmark or external control cohort has complete ascertainment of death.

## Summary

Single-arm phase II trials in glioblastoma are interpreted against historical benchmarks or external controls drawn from earlier trials or real-world data. Trial patients who stop treatment because of AEs may leave survival follow-up, whereas deaths in registry-based or institutional cohorts are ascertained regardless of treatment. When AE-related discontinuation identifies patients at higher risk of death, censoring at discontinuation removes high-risk patients from the trial cohort only.

The simulation study follows the ADEMP structure (Morris, White and Crowther, Stat Med 2019): 60 scenarios with 5,000 replicates each, 60 trial patients compared with 180 external controls in recurrent glioblastoma (control median overall survival 9 months), varying the cumulative incidence of AE-related discontinuation (10% to 30%), the hazard ratio for death after discontinuation (1 to 3), the probability of loss to follow-up at discontinuation (0.5 or 1), the loss mechanism, the sample size and the disease setting. Estimands follow the treatment-policy strategy. Naive analyses, inverse probability of censoring weighting on baseline covariates or on discontinuation status, and delta-adjusted multiple imputation (tipping point) are compared on the hazard-ratio scale and on clinical scales (median overall survival, 12-month overall survival, benchmark comparison).

Main findings for an ineffective treatment with 20% discontinuation, a doubled post-discontinuation hazard and follow-up ending at discontinuation:

- naive hazard ratio 0.879 with a two-sided type I error of 10.2%;
- median overall survival overestimated by 1.2 months and 12-month overall survival by 5.5 percentage points;
- benchmark comparison of 12-month overall survival positive in 11.2% of trials (nominal 2.5%);
- weighting on baseline covariates removed at most 18% of the bias; weighting on discontinuation status removed it when some discontinuing patients remained under follow-up.

A review of 55 single-arm phase II or I/II trials in high-grade glioma (PubMed Central, 2015 to 2025) found overall survival as the primary endpoint in 19, interpretation against historical benchmarks in 22, numerical reporting of AE-related discontinuation in 17, a statement of continued survival follow-up after discontinuation in 3 and an informative-censoring sensitivity analysis in 0.

## Repository structure

```
R/functions.R                   data-generating mechanism, calibration, estimators, performance measures
R/scenarios.R                   scenario definitions
R/01_calibrate_truth.R          calibration and true values (10^6 patients per cohort)
R/02_run_simulation.R           Monte Carlo simulation (scenario-specific seeds)
R/03_tipping_point.R            delta-adjusted multiple imputation
R/04_summarise.R                performance measures with Monte Carlo standard errors
R/05_figures.R                  figure code
review/01_search_screen.py      PubMed Central search and full-text screening
review/02_summarise_review.py   summary of the coded review
review/03_pubmed_denominator.py abstract-level PubMed query for representativeness
review/review_coded.csv         manually verified coding of 63 screened articles (55 eligible)
review/review_summary.csv
review/pubmed_denominator.csv
results/scenarios_truth.csv     calibrated parameters and true values
results/performance.csv         performance of all estimators in all scenarios
results/tipping_point_performance.csv
run_all.sh
session_info.txt
```

## Reproduction

Requirements: R 4.3 or later with `survival`, `ggplot2` and `patchwork` (and base `parallel`); Python 3 with `pandas` for the review scripts. The review scripts use the NCBI E-utilities; an API key can be supplied in the environment variable `NCBI_API_KEY`.

```
./run_all.sh
```

The full pipeline takes approximately 75 minutes on two cores. Replicate-level estimates are written to `results/raw/` (about 30 MB; not tracked) and are reproduced exactly from the scenario-specific seeds, independent of the number of cores. `R/02_run_simulation.R` skips scenarios whose output already exists.

## Licence

MIT. See `LICENSE`.

## Citation

See `CITATION.cff`.
