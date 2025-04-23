#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=flare_model
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16gb
#SBATCH --output=logs/summary.%j.log
#SBATCH --time=24:00:00

set -e

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility

module purge
module list

## Job commands here
CHROM=1
THREADS=${SLURM_CPUS_PER_TASK}
OUTPUT_PREFIX="LIBD_TOPMed_AA_EA"
SOFTWARE="/projects/p32505/opt/bin"
BASE_DIR="/projects/b1213/resources"
MAP_DIR="${BASE_DIR}/1kGP/genetic_maps"
REFERENCE_DIR="${BASE_DIR}/1kGP/GRCh38_phased_vcf"
QUERY_DIR="${BASE_DIR}/libd_data/genotypes/combined_data/AA_EA/flare-vcf-no-missing"

log_message "**** FLARE Local Ancestry Analysis ****"
# Check if required files exist
if [ ! -d "$REFERENCE_DIR" ]; then
    log_message "Error: Reference directory not found: $REFERENCE_DIR"
    exit 1
fi

if [ ! -d "$QUERY_DIR" ]; then
    log_message "Error: Query directory not found: $QUERY_DIR"
    exit 1
fi

REF_VCF="1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
QUERY_VCF="TOPMed_LIBD.AA_EA.chr${CHROM}.vcf.gz"

java -Xmx16g -jar $SOFTWARE/flare.jar \
     ref="${REFERENCE_DIR}/raw/${REF_VCF}" \
     ref-panel="${REFERENCE_DIR}/local-ancestry-ref/samples_id2" \
     map="${MAP_DIR}/plink.chr${CHROM}.GRCh38.map" \
     gt="${QUERY_DIR}/_m/${QUERY_VCF}" nthreads=$THREADS \
     seed=13131313 array=true out="${OUTPUT_PREFIX}"

log_message "**** Job ends ****"
