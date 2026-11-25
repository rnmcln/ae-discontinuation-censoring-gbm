#!/usr/bin/env bash
# Full pipeline. Approximate run time on two cores: 75 minutes.
set -euo pipefail
Rscript R/01_calibrate_truth.R
Rscript R/02_run_simulation.R 5000 2
Rscript R/03_tipping_point.R 1000 2
Rscript R/04_summarise.R
LC_ALL=C.UTF-8 Rscript R/05_figures.R figures
( cd review && python3 02_summarise_review.py )
