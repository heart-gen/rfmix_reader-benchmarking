#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=index_chr
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.jr1@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22
#SBATCH --mem=20gb
#SBATCH --output=log_files/indexing.%A_%a.log
#SBATCH --time=00:30:00

# Function to echo with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"


log_message "**** Quest info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

## List current modules for reproducibility
log_message "**** Loading modules ****"

module purge
module load htslib/1.16
module list

CHROM=${SLURM_ARRAY_TASK_ID}
VCF_FILE="chr${CHROM}.vcf.gz"
TEMP_FILE="chr${CHROM}.vcf.gz.tmp"

log_message "Run job for ground-truth chromosome: ${CHROM}"

# Function to check BGZF
is_bgzf() {
    # Returns 0 if BGZF, 1 if not
    if [[ $(htsfile "$1" 2>/dev/null) == *"BGZF"* ]]; then
        return 0
    else
        return 1
    fi
}

log_message "*** Checking compression format for $VCF_FILE ***"
if ! is_bgzf "$VCF_FILE"; then
    log_message "$VCF_FILE is not BGZF — re-compressing..."
    gunzip -c "$VCF_FILE" | bgzip -c > "$TEMP_FILE" && mv "$TEMP_FILE" "$VCF_FILE"
else
    log_message "$VCF_FILE is already BGZF — skipping recompression."
fi

log_message "*** Indexing VCF ***"
tabix -p vcf "$VCF_FILE"

log_message "**** Job ends ****"
