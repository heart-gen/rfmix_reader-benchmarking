#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=convert_references
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=32
#SBATCH --array=1-22
#SBATCH --time=12:00:00
#SBATCH --output=logs/reference_conversion.%A_%a.log

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

OUTDIR="reference_zarr"
CHR=${SLURM_ARRAY_TASK_ID}
CHR_OUT="${OUTDIR}/chr${CHR}.zarr"
VCF_DIR="/ocean/projects/bio250020p/shared/resources/1kGP/GRCh38_phased_vcf/raw"
VCF="${VCF_DIR}/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

CHUNK_LEN=100000
SAMPLE_CHUNK=2048
WORKERS=32

mkdir -p "${OUTDIR}"

log_message "Processing chromosome ${CHR}"
log_message "VCF: ${VCF}"
log_message "Zarr output: ${CHR_OUT}"

# Run prepare-reference
prepare-reference \
    --chunk-length ${CHUNK_LEN} \
    --samples-chunk-size ${SAMPLE_CHUNK} \
    --worker-processes ${WORKERS} \
    --verbose "${CHR_OUT}" "${VCF}"

conda deactivate
log_message "Job finished at: $(date)"
