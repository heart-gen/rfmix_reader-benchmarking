#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=gen_plinks
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=32
#SBATCH --time=01:00:00
#SBATCH --output=logs/plink_merge.%J.log

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
OUTDIR="plink-files"
mkdir -p ${OUTDIR}

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/genomics

ls plink-files/chr*.pgen | sort -V | sed 's/.pgen//' > merge_list.txt

plink2 --pmerge-list merge_list.txt --make-pgen --out "${OUTDIR}/simulated"

conda deactivate
log_message "**** Job ends ****"
