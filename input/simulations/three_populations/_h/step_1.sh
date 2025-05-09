#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --job-name=simu_3pops
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=150gb
#SBATCH --output=logs/output.%A_%a.log
#SBATCH --time=03:00:00

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
module load htslib/1.16
module list

## Edit with your job command
log_message "**** Loading mamba environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/env/AI_env

echo "**** Run simulation ****"
ONE_K="/projects/b1213/resources/1kGP/"
CHROM=${SLURM_ARRAY_TASK_ID}
INPUTS="simulation-files"
VCFDIR="vcf-files"

haptools simgenotype \
         --model ${INPUTS}/AFR_washington.dat \
         --mapdir ${ONE_K}/genetic_maps/ \
         --chroms ${CHROM} \
         --seed 20240126 \
         --ref_vcf ${ONE_K}/GRCh38_phased_vcf/raw/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
         --sample_info ${INPUTS}/1k_sampleinfo.tsv \
         --out ${VCFDIR}/chr${CHROM}.vcf.gz

tabix -p vcf ${VCFDIR}/chr${CHROM}.vcf.gz
mv ${VCFDIR}/*bp ${INPUTS}/

conda deactivate
log_message "**** Job ends ****"
