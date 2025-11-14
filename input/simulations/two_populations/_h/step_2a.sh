#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --job-name=rfmix_simu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=250gb
#SBATCH --output=logs/rfmix.%A_%a.log
#SBATCH --array=1-22
#SBATCH --time=128:00:00

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
OUTDIR="rfmix-out"
VCFDIR="simulation-files"
CHROM=${SLURM_ARRAY_TASK_ID}
SOFTWARE="/projects/p32505/opt/bin"
ONE_K="/projects/b1213/resources/1kGP/"
REF="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/local-ancestry-ref"

log_message "**** Run RFMix ****"
echo -e "Chromosome: ${CRHOM}"

mkdir -p $OUTDIR

$SOFTWARE/rfmix \
    -f ${VCFDIR}/chr${CHROM}.vcf.gz \
    -r $REF/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz \
    -m $REF/samples_id2 \
    -g $REF/genetic_map38 \
    -o ${OUTDIR}/chr${CHROM} \
    --chromosome=chr${CHROM}

log_message "**** Job ends ****"
