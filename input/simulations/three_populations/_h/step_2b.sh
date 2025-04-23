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
#SBATCH --output=prep.%A_%a.log
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
module list

## Edit with your job command
CHROM=${SLURM_ARRAY_TASK_ID}
VCFDIR="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/raw/"

log_message "**** Prepare samples ****"
echo -e "Chromosome: ${CRHOM}"

bcftools view -v snps -S ./temp/samples_id -Oz \
         -o ./temp/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz \
         $VCFDIR/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz

log_message "**** Job ends ****"
