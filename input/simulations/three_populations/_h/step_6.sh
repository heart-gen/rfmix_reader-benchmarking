#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=patch_contigs
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=25gb
#SBATCH --output=logs/patch.%A_%a.log
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
module load bcftools/1.10.1
module list

## Edit with your job command
echo "**** Run simulation ****"
GTDIR="gt-files"
CHROM=${SLURM_ARRAY_TASK_ID}
TEMPFILE="temp/temp.chr${CHROM}.vcf.gz"
REF="/projects/b1213/resources/1kGP/references/Homo_sapiens_assembly38.fasta.fai"

# Run haptools to generate ground truth data
VCF="${GTDIR}/chr${CHROM}.vcf.gz"
mv -v ${VCF} "${TEMPFILE}"

bcftools annotate --contigs "$REF" -O z -o "${VCF}" "${TEMPFILE}"

# Index output
rm ${VCF}.tbi
tabix -p vcf ${VCF}

rm $TEMPFILE

log_message "**** Job ends ****"
