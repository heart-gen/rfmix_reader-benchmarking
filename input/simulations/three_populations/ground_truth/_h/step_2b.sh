#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=convert_int
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.jr1@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40gb
#SBATCH --output=log_files/convert.%A_%a.log
#SBATCH --array=4
#SBATCH --time=4:00:00

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
module load gcc/12.3.0-gcc
module list

## Edit with your job command
log_message "**** Loading mamba environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/env/AI_env

log_message "**** Run conversion ****"
CHROM=${SLURM_ARRAY_TASK_ID}

python ../_h/02.convert_hap_to_geno.py ${CHROM} \
       per_snp_ancestry.chr${CHROM}.tsv.gz chr${CHROM}.vcf.gz 10000

conda deactivate
log_message "**** Job ends ****"
