#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=fix_files
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --array=1-22
#SBATCH --output=logs/fix.%A_%a.log
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

module purge
module list

## Edit with your job command

# path to your contigs file
CONTIGS="../../three_populations/_m/contigs.txt"
VCFDIR="gt-files"
CHR="chr${SLURM_ARRAY_TASK_ID}"
OUT="${VCFDIR}/${CHR}.vcf.gz"
IN="${VCFDIR}/back/${CHR}.vcf.gz"

log_message "Processing ${CHR} ..."

# grab correct contig line from contigs.txt
CONTIG_LINE=$(grep -w "ID=${CHR}" "$CONTIGS")
if [[ -z "$CONTIG_LINE" ]]; then
    echo "  ERROR: No contig line found for ${CHR} in $CONTIGS"
    continue
fi

log_message "**** Loading conda environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/genomics

# extract header, replace old contig line with corrected one
bcftools view -h "$IN" \
    | sed "s/^##contig=<ID=${CHR}>.*/${CONTIG_LINE}/" > header.${CHR}.tmp

# reheader the file
bcftools reheader -h header.${CHR}.tmp -o "$OUT" "$IN"

# index
tabix -p vcf "$OUT"

echo "  Fixed header and reindexed: $OUT"
log_message "**** Job ends ****"

