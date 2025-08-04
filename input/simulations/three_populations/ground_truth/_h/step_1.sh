#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=convert_bp
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --output=log_files/output.%A_%a.log
#SBATCH --array=1-22
#SBATCH --time=01:00:00

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

log_message "**** Run conversion ****"
INPUT="../../_m/simulation-files"
CHROM=${SLURM_ARRAY_TASK_ID}

python ../_h/01.convert_bp_to_snp.py \
       ${INPUT}/chr${CHROM}.bp \
       ${INPUT}/chr${CHROM}.vcf.gz \
       per_snp_ancestry.chr${CHROM}.tsv.gz

conda deactivate
log_message "**** Job ends ****"
