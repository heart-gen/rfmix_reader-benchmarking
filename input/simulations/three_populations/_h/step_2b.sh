#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=prep_data
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=24
#SBATCH --array=1-22
#SBATCH --output=logs/prep.%A_%a.log
#SBATCH --time=01:30:00

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module load anaconda3/2024.10-1
module list

## Edit with your job command
CHROM=${SLURM_ARRAY_TASK_ID}
TEMPDIR="temp"
FILTERED_VCF="${TEMPDIR}/chr${CHROM}.biallelic.vcf.gz"

log_message "**** Prepare samples ****"
echo -e "Chromosome: ${CHROM}"

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics

bcftools view -v snps -S "${TEMPDIR}/samples_id" -Oz \
         -o "${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz" \
         "${FILTERED_VCF}"

tabix -f -p vcf "${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz"

conda deactivate
log_message "**** Job ends ****"
