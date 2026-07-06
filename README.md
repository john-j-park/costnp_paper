# Flexible Neyman–Pearson Classification via Cost-Sensitive Learning

**John Park, Rachel Wang, Jessica Li, Xin Tong**

Code to reproduce all experiments and figures in the paper.

## Dependencies

All code is written in R. The following packages are required:

```r
install.packages(c("ggplot2", "dplyr", "data.table", "furrr", "future", "MASS", "mvtnorm", "e1071"))
```

The `npcost` and `nproc` packages are also required:

```r
# nproc is available on CRAN
install.packages("nproc")

# npcost — install from source (link TBD)
```

## Reproducing the paper

All scripts should be run from the **project root** (the directory containing `scripts/` and `data/`). Each experiment writes results to `output/`; move the chosen run into `data/` before running the corresponding plot script (see workflow below).

### Figures in Section 2

```r
source("scripts/plot_section_2.R")
```

### Introduction figure

```r
source("scripts/plot_introduction.R")
```

### Experiments

Run each experiment script, then its corresponding plot script.

| Experiment | Run | Plot |
|---|---|---|
| Experiment 1 | `experiment_1.R` | `plot_experiment_1_main.R`, `plot_experiment_1_supplement.R` |
| Experiment 2 | `experiment_2.R` | `plot_experiment_2.R` |
| Experiment 3 | `experiment_3.R` | `plot_experiment_3.R` |
| Experiment 4 | `experiment_4.R` | `plot_experiment_4.R` |
| Real data    | `experiment_real_data.R` | `plot_experiment_real_data.R`, `latex_real_data.R` |

**Workflow for each experiment:**
1. Run `scripts/experiment_X.R` — results are saved to `output/` with a timestamp.
2. Move and rename the output file into `data/` (e.g. `data/experiment_1.RData`).
3. Run the corresponding plot script.

Pre-computed results for all experiments are included in `data/` so figures can be reproduced without re-running the simulations.

### Shared utilities

`scripts/paper_utils.R` is sourced automatically by the other scripts — it does not need to be run directly.

## Data

The `data/` folder contains both raw datasets and pre-computed experiment results.

| File | Source |
|---|---|
| `CTG.csv` | [UCI Cardiotocography dataset](https://archive.ics.uci.edu/dataset/193/cardiotocography) |
| `diabetes.csv` | [UCI Diabetes dataset](https://archive.ics.uci.edu/dataset/34/diabetes) |
| `mammography.csv` | [OpenML Mammography dataset](https://www.openml.org/d/310) |
| `phoneme.csv` | [OpenML Phoneme dataset](https://www.openml.org/d/1489) |
| `winequality-white.csv` | [UCI Wine Quality dataset](https://archive.ics.uci.edu/dataset/186/wine+quality) |
| `magic04.data` | [UCI MAGIC Gamma Telescope dataset](https://archive.ics.uci.edu/dataset/159/magic+gamma+telescope) |
| `transfusion.data` | [UCI Blood Transfusion dataset](https://archive.ics.uci.edu/dataset/176/blood+transfusion+service+center) |
| `*.RData` | Pre-computed experiment results |
