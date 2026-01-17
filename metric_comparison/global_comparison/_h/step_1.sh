#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=global_ancestry
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --output=global-ancestry.%j.log
#SBATCH --time=01:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility
module list

echo "**** Run global ancestry analysis ****"

SIM_INPUT=${SIM_INPUT:-"/path/to/simulation/input"}
RFMIX_INPUT=${RFMIX_INPUT:-"/path/to/rfmix/input"}
FLARE_INPUT=${FLARE_INPUT:-"/path/to/flare/input"}
OUTPUT_DIR=${OUTPUT_DIR:-"./outputs/global_ancestry"}

python ../_h/per_chrom_global_ancestry.py \
    --simu-input "${SIM_INPUT}" \
    --rfmix-input "${RFMIX_INPUT}" \
    --flare-input "${FLARE_INPUT}" \
    --output "${OUTPUT_DIR}"

echo "**** Job ends ****"
date
