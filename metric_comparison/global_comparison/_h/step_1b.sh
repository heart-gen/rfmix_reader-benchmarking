#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=two_global
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=40
#SBATCH --time=08:00:00
#SBATCH --output=logs/global-ancestry.two_pop.%j.log

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
OUTPUT_DIR="results"
SIMU_DIR="input/simulations/two_populations/_m/gt-files"
RFMIX_DIR="input/simulations/two_populations/_m/rfmix-out"
FLARE_INPUT="input/simulations/two_populations/_m/flare-out"

python ../_h/01.global_ancestry.py \
    --simu-input "${SIM_INPUT}" --rfmix-input "${RFMIX_INPUT}" \
    --flare-input "${FLARE_INPUT}" --output "${OUTPUT_DIR}" --population "two"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
