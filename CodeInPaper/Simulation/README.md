# Simulation: calling density-correlated genes of known status

Reproduces the calling-rule simulation of the Results subsection "Validation of
DenMark on simulated data" (Figure S3; Supplemental Information, "Simulation:
Identification of Density-Correlated Genes"). Each gene is simulated under the
DenMark model with a known correlation rho between the cell-density and
expression fields, fitted with `DenMark.stan`, and called a density-correlated
gene (DCG) when the 95% credible interval of rho (`Astar[2,1]`) excludes zero.

## Scripts

| File | Purpose | Output (below `OUTPUT_DIR`) |
|---|---|---|
| `00_functions.R` | Simulates cell counts and expression on a grid from two Matérn 3/2 fields with correlation rho, assembles the Stan data list, fits the model and applies the calling rule. Sourced by the two experiments. | none |
| `02_experiment_FPR.R` | Null genes (rho = 0): false positive rate | `Revision-R1/FPR+TPR/fpr_results.rds`, `fpr_histogram.png` (Figure S3B) |
| `04_experiment_TPR.R` | Non-null genes (rho = 0.5 and -0.5, half each): true positive rate | `Revision-R1/FPR+FDR/tpr_results.rds`, `tpr_histogram.png` (Figure S3A) |

## Inputs

None: the data are simulated. The model is `inst/stan/DenMark.stan` of the
checkout, located through `REPO_ROOT`.

## How to run

Set `OUTPUT_DIR` and `REPO_ROOT` in the CONFIGURE ME block of each experiment and
run it from the repository root, for example
`Rscript CodeInPaper/Simulation/02_experiment_FPR.R`. The two experiments are
independent. The number of genes (`N_SIM`), the grid (`make_grid()`) and the
number of basis functions per dimension (`M_BASIS`) are set near the top of each
experiment; the Supplemental Information gives the values behind Figure S3.

## Compute notes

Each gene is one fit of 2 chains with 400 warm-up and 400 sampling iterations.
The fits run in parallel through `parallel::mclapply()`, which forks on Linux and
macOS and runs serially on Windows (the FPR script describes a `future.apply`
alternative). Fits with R-hat of 1.05 or more are dropped before the rates are
computed.
