#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=gen_plinks
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=32
#SBATCH --array=1-22
#SBATCH --time=01:00:00
#SBATCH --output=logs/plink.%A_%a.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility

module purge
module load anaconda3/2024.10-1
module list

## Edit with your job command

# path to your contigs file
VCFDIR="gt-files"
OUTDIR="plink-files"
CHR=${SLURM_ARRAY_TASK_ID}

mkdir -p ${OUTDIR}

log_message "Processing ${CHR} ..."

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics

plink2 --vcf "${VCFDIR}/chr${CHR}.vcf.gz" --make-pgen --out "${OUTDIR}/chr${CHR}"

conda deactivate
log_message "**** Job ends ****"
