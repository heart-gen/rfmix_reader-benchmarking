#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=impute_three_pop
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=64
#SBATCH --time=12:00:00
#SBATCH --array=0-2
#SBATCH --output=logs/impute.three_pop.%A_%a.log

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

METHODS=("linear" "nearest" "stepwise")
METHOD="${METHODS[${SLURM_ARRAY_TASK_ID}]}"

if [ -z "${METHOD}" ]; then
    echo "Invalid method index: ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/ai_env

log_message "**** Run analysis ****"
RFMIX_DIR="input/simulations/three_populations/_m/rfmix-files"

python ../_h/01.impute_data.py \
       --rfmix-input "${RFMIX_DIR}" \
       --population "three" --method "${METHOD}"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
