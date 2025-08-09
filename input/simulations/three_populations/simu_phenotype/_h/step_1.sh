#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --job-name=simu_pheno
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
module list

## Edit with your job command
log_message "**** Loading mamba environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/env/AI_env

echo "**** Run simulation ****"
CHROM=${SLURM_ARRAY_TASK_ID}
OUTDIR="haplotype-files"
VCFDIR="../../_m/vcf-files"
PGEN="plink-files"

mkdir -p ${OUTDIR}
mkdir -p ${PGEN}

haptools transform \
         --ancestry --verbosity INFO \
         --output ${PGEN}/chr${CHROM}.pgen \
         ${VCFDIR}/chr${CHROM}.vcf.gz \
         ${OUTDIR}/chr${CHROM}.hap

conda deactivate
log_message "**** Job ends ****"
