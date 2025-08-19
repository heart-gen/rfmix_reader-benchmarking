#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=viz_local
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.jr1@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --output=logs/global_ancestry.%J.log
#SBATCH --time=00:30:00

# Function to echo with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module load gcc/12.3.0-gcc
module list

## Edit with your job command
log_message "**** Loading mamba environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/env/AI_env

python ../_h/05.global_ancestry-ground_truth.py \
  --folder /projects/p32505/users/manuel/rfmix_reader-benchmarking/input/simulations/two_populations/ground_truth/_m/ \
  --file global_ancestry.tsv \
  --chromosome_plots 

if [ $? -ne 0 ]; then
    log_message "Error: mamba or script execution failed"
    exit 1
fi

conda deactivate
log_message "**** Job ends ****"
