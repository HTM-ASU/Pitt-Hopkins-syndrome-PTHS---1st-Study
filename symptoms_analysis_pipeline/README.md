# FLACC / DSR reproducible analysis pipeline

This repository reproduces the analysis for three outcomes:
- FLACC
- DSR Total Events
- GSRS AVG

It uses a single input Excel file and a single configuration file so the code does not need to be edited when the data file stays the same.

## Required input file
Place the Excel file here:

`data/raw/Final_merged_data_for_github_test_dataset.xlsx`

Use this sheet name:

`Merged`

Required columns in the input file:

- `Participant`
- `Group`
- `Timepoint`
- `FLACC`
- `DSR Total Events`
- `GSRS AVG`

Expected timepoints:

- `Baseline`
- `End Part 1`
- `End Part 2`
- `End Part 3`

Expected groups:

- `GROUP A`
- `GROUP B`

## What the pipeline does
1. Computes per-patient change scores as `Baseline - follow-up`.
2. Runs Welch one-sided t-tests comparing `GROUP A` vs `GROUP B` for each score and each follow-up timepoint.
3. Creates boxplots with jittered points for all three scores.
4. Adds only the `End Part 1` p-value label above the corresponding panel.

## Files to run
- `run_all.sh`
- `scripts/python/01_prepare_change_scores.py`
- `scripts/python/02_run_ttests.py`
- `scripts/R/03_plot_results.R`

## Setup
### Python
Create a clean environment. Recommended:

```bash
conda env create -f environment.yml
conda activate flacc_dsr
```

If you prefer pip, use:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### R packages
Install once in R:

```r
install.packages(c("readxl", "yaml", "tidyr", "ggplot2", "dplyr"))
```

## Run everything
From the repository root:

```bash
bash run_all.sh
```

## Outputs
After the run finishes, look in:

- `output/tables/Per_patient_differences_FLACC_DSR.xlsx`
- `output/tables/T-test_results_FLACC_DSR_greater_in_treated.xlsx`
- `output/figures/Symptoms_FLACC_DSR_Delta_p-values_04.06.2026_bold.pdf`
- `output/logs/`

## Notes
- The plotting script keeps `Baseline` on the x-axis.
- No connecting lines are drawn for the p-value annotation.
- The p-value label is shown only for `End Part 1`.
