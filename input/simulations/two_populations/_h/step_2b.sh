#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=gzip_rfmix
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --output=logs/gzip_rfmix.%A_%a.log
#SBATCH --time=00:30:00

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
CHROM=${SLURM_ARRAY_TASK_ID}
OUTDIR="rfmix-files"

log_message "**** Run compression ****"
echo -e "Chromosome: ${CRHOM}"

# Files to compress
for suffix in fb msp sis; do
    INFILE=${OUTDIR}/chr${CHROM}.${suffix}.tsv
    if [[ -f $INFILE ]]; then
        echo "Compressing $INFILE"
        gzip -f "$INFILE"
    else
        echo "Warning: $INFILE not found"
    fi
done
log_message "**** Job ends ****"
