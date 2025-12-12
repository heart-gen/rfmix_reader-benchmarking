#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=phase_merged_data
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=64
#SBATCH --array=0-1
#SBATCH --time=12:00:00
#SBATCH --output=logs/merge_phased.%A_%a.log

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
conda activate /ocean/projects/bio250020p/shared/opt/env/ai_env

log_message "**** Run analysis ****"
RFMIX_LABELS=("RFMIX_TWO" "RFMIX_THREE")
RFMIX_PATHS=(
    "../../input/simulations/two_populations/_m/rfmix-out/phased_files"
    "../../input/simulations/three_populations/_m/rfmix-files/phased_files"
)

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
PHASED_DIR=${RFMIX_PATHS[$TASK_ID]}
DATASET_LABEL=${RFMIX_LABELS[$TASK_ID]}

if [ -z "${PHASED_DIR}" ]; then
    echo "Invalid or missing SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

log_message "Processing dataset: ${DATASET_LABEL} (${PHASED_DIR})"

declare -a CHR_FILES
for chr in {1..22}; do
    CHR_FILES+=("${PHASED_DIR}/phased_chr${chr}.zarr")
done

merge-phased-zarrs "${PHASED_DIR}/phased_data.zarr" "${CHR_FILES[@]}"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
