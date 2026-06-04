#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=rfmix_msp_cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=16
#SBATCH --array=1-3
#SBATCH --time=04:00:00
#SBATCH --output=logs/rfmix_msp_cpu.%A_%a.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

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
conda activate rfmix-reader-msp-cpu

TASK="${SLURM_ARRAY_TASK_ID}"
INPUT_DIR="${RFMIX_MSP_DIR:-../../../input/aanri_data/rfmix-version/_m}"
OUTDIR="aanri_msp_cpu"

python ../_h/01.msp_target_query.py --input "${INPUT_DIR}" \
       --output "${OUTDIR}" --label "task_${TASK}" --task "${TASK}" \
       --target-count "${RFMIX_MSP_TARGET_COUNT:-1000}"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
