#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=prep_data
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --output=logs/prep.%A_%a.log
#SBATCH --time=01:30:00

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
log_message "**** Loading modules ****"

module purge
module load bcftools/1.10.1
module load htslib/1.16
module list

## Edit with your job command
CHROM=${SLURM_ARRAY_TASK_ID}
TEMPDIR="temp"
FILTERED_VCF="${TEMPDIR}/chr${CHROM}.biallelic.vcf.gz"

log_message "**** Prepare samples ****"
echo -e "Chromosome: ${CHROM}"

bcftools view -v snps -S "${TEMPDIR}/samples_id" -Oz \
         -o "${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz" \
         "${FILTERED_VCF}"

tabix -f -p vcf "${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz"
log_message "**** Job ends ****"
