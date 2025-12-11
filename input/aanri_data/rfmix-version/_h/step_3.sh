#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=phase_data
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=32
#SBATCH --time=4:00:00
#SBATCH --array=1-22
#SBATCH --output=logs/phase.annri_data.%A_%a.log

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
CHR=${SLURM_ARRAY_TASK_ID}
RFMIX_DIR="input/aanri_data/rfmix-version/_m"
REF_DIR="input/references/_m/two_populations/reference_zarr/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.zarr"

python ../_h/03.phase_data.py \
       --rfmix-input "$RFMIX_DIR" \
       --chrom ${CHR} --ref-input "$REF_DIR"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
