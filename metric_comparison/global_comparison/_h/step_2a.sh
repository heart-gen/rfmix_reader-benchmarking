#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=phased_global_three_pop
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=16
#SBATCH --time=03:00:00
#SBATCH --array=1-22
#SBATCH --output=logs/phased_global.three_pop.%A_%a.log

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
CHR=${SLURM_ARRAY_TASK_ID}
OUTPUT_DIR="phased"
SIMU_DIR="input/simulations/three_populations/_m/gt-files"
RFMIX_DIR="input/simulations/three_populations/_m/rfmix-files"
FLARE_DIR="input/simulations/three_populations/_m/flare-out"
SAMPLE_ANNOT="input/references/_m/three_populations/reference_zarr/samples_id2"
REF_DIR="input/references/_m/three_populations/reference_zarr"

python ../_h/02.phased_global_ancestry.py \
       --rfmix-input "$RFMIX_DIR" \
       --simu-input "$SIMU_DIR" \
       --flare-input "$FLARE_DIR" \
       --sample-annot "$SAMPLE_ANNOT" \
       --ref-input "$REF_DIR" \
       --output "$OUTPUT_DIR" \
       --population "three" \
       --chrom "$CHR"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
