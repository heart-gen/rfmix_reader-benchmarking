#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=flare_model
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24gb
#SBATCH --output=logs/summary.%j.log
#SBATCH --time=02:00:00

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
SOFTWARE="/projects/p32505/opt/bin"
BASE_DIR="/projects/b1213/resources"
MAP_BASE="/projects/p32505/projects/rfmix_reader-benchmarking"
MAP_DIR="${MAP_BASE}/input/simulations/three_populations/_m/temp"
REFERENCE_DIR="${BASE_DIR}/1kGP/GRCh38_phased_vcf"
QUERY_DIR="${BASE_DIR}/libd_data/genotypes/combined_data/allsamples/vcf-format"

log_message "**** FLARE Local Ancestry Analysis ****"

REF_VCF="1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
QUERY_VCF="chr${CHROM}_phased.vcf.gz"

java -Xmx24g -jar $SOFTWARE/flare.jar \
     ref="${REFERENCE_DIR}/raw/${REF_VCF}" \
     ref-panel="${REFERENCE_DIR}/local-ancestry-ref/samples_id2" \
     map="${MAP_DIR}/plink.chr${CHROM}.GRCh38.reformatted.map" \
     gt="${QUERY_DIR}/_m/${QUERY_VCF}" nthreads=$THREADS \
     seed=13131313 array=true out="chr${CHROM}"

log_message "**** Job ends ****"
