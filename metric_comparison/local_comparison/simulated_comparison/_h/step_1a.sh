#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=unphased_metrics_three_pop
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=24
#SBATCH --time=08:00:00
#SBATCH --array=5-22
#SBATCH --output=logs/unphased_metrics.three_pop.%A_%a.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/ml_dev

log_message "**** Run analysis ****"
SIMU_DIR="input/simulations/three_populations/_m/gt-files"
RFMIX_DIR="input/simulations/three_populations/_m/rfmix-files"
CHR=${SLURM_ARRAY_TASK_ID}

python ../_h/01.unphased_simulation.py \
       --rfmix-input "$RFMIX_DIR" --simu-input "$SIMU_DIR" \
       --output "unphased" --population "three" --chrom "$CHR"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
