#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=convert_references
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=64
#SBATCH --array=2
#SBATCH --time=24:00:00
#SBATCH --output=logs/reference_conversion.two_pop.%A_%a.log

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

CHR=${SLURM_ARRAY_TASK_ID}
OUTDIR="two_populations/reference_zarr"
VCF_DIR="/ocean/projects/bio250020p/shared/resources/1kGP/GRCh38_phased_vcf"
VCF="${VCF_DIR}/raw/1kGP_high_coverage_Illumina.chr${CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"

WORKERS=32
CHUNK_LEN=250000
SAMPLE_CHUNK=2048

mkdir -p "${OUTDIR}"

log_message "Processing chromosome ${CHR}"
log_message "VCF: ${VCF}"
log_message "Zarr output: ${OUTDIR}"

# Run prepare-reference
prepare-reference \
    --chunk-length ${CHUNK_LEN} \
    --samples-chunk-size ${SAMPLE_CHUNK} \
    --worker-processes ${WORKERS} \
    --verbose "${OUTDIR}" "${VCF}"

cp -v "${VCF_DIR}/local-ancestry-ref/samples_id2" "${OUTDIR}/samples_id2"

conda deactivate
log_message "Job finished at: $(date)"
