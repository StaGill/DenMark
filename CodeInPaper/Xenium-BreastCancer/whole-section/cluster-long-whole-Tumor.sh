#!/bin/bash
#SBATCH --array=1-313            # upper bound; tasks past the QC-passed count exit at once
#SBATCH --account=def-<account>  # replace this with your supervisors account
#SBATCH --ntasks=10             # number of processes
#SBATCH --mem-per-cpu=4000M     # memory; default unit is megabytes
#SBATCH --time=05:00:00         # time (HH:MM:SS)
#SBATCH --output=logs/%x-%j.out

# ---- CONFIGURE ME -------------------------------------------------------------
# Exported for DenMark-BC-RA-cluster-whole.R, which reads them with Sys.getenv().
# Relative paths resolve from the submission folder (this folder). Create logs/
# here before submitting (mkdir -p logs), since #SBATCH lines cannot use variables.
export DATA_DIR=../../../data      # holds Dataset/xenium-breast-cancer/rep1/
export OUTPUT_DIR=../../../output  # holds qc-whole/ from qc-select-genes-whole.R; results go to BC-whole/
export REPO_ROOT=../../..          # root of the DenMark checkout
# ------------------------------------------------------------------------------

# Modules
module load StdEnv/2020 gcc/9.3.0 r/4.3.1

# Export the node names.
export NODESLIST=$(echo $(srun hostname | cut -f 1 -d '.'))

# Group 4 = Tumor (DCIS_1 + DCIS_2 + Invasive_Tumor)   (division C)
GROUP=4

R -f DenMark-BC-RA-cluster-whole.R --args $SLURM_ARRAY_TASK_ID $GROUP

# Efficient submit:
#   sbatch --array=1-$(cat $OUTPUT_DIR/qc-whole/Tumor_ngenes.txt) cluster-long-whole-Tumor.sh
