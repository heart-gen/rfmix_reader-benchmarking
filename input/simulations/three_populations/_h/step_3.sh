#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=long
#SBATCH --job-name=rfmix_3pop
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=175gb
#SBATCH --output=logs/rfmix.%A_%a.log
#SBATCH --time=72:00:00

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
SOFTWARE="/projects/p32505/opt/bin"
OUTDIR="rfmix-files"
VCFDIR="vcf-files"
INPUTS="temp"

log_message "**** Run RFMix ****"
echo -e "Chromosome: ${CRHOM}"

$SOFTWARE/rfmix \
    -f ${VCFDIR}/chr${CHROM}.vcf.gz \
    -r ${INPUTS}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz \
    -m ${INPUTS}/samples_id2 \
    -g ${INPUTS}/genetic_map38 \
    -o ${OUTDIR}/chr${CHROM} \
    --chromosome=chr${CHROM}

log_message "**** Job ends ****"
