#!/bin/bash
#SBATCH --account=def-<account>     # CONFIGURE ME: your allocation
#SBATCH --cpus-per-task=4           # one core per chain (chains = parallel_chains = 4)
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=05:00:00             # per gene; the paper's whole-section Xenium scan asked for 5 hours
#SBATCH --output=logs/%x-%A_%a.out
#
# One DenMark gene per SLURM array task. Copy this file and denmark-scan.R into a
# working folder (they are in system.file("slurm", package = "DenMark")), set the
# CONFIGURE ME lines in both, then run from that folder:
#
#   mkdir -p logs
#   Rscript denmark-scan.R prepare
#   sbatch --array=1-$(cat denmark_scan/n_genes.txt) denmark-scan-array.sh
#   Rscript denmark-scan.R collect      # after the array has finished

module load r                       # CONFIGURE ME: an R with DenMark, cmdstanr and CmdStan installed

Rscript denmark-scan.R gene "$SLURM_ARRAY_TASK_ID"
