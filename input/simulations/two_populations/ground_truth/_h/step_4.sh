#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=long
#SBATCH --job-name=rfmix_ground
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.jr1@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb
#SBATCH --output=log_files/rfmix.%A_%a.log
#SBATCH --array=1-22
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
OUTDIR="rfmix-out"
SOFTWARE="/projects/p32505/opt/bin"
ONE_K="/projects/b1213/resources/1kGP/"
REF="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/local-ancestry-ref"

mkdir -p $OUTDIR

log_message "**** Run RFMix ****"
echo -e "Chromosome: ${CHROM}"

# Set thread count for parallelization
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
THREADS=${SLURM_CPUS_PER_TASK}

$SOFTWARE/rfmix \
    -f chr${CHROM}.vcf.gz \
    -r $REF/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz \
    -m $REF/samples_id2 \
    -g $REF/genetic_map38 \
    -o ${OUTDIR}/chr${CHROM} \
    --n-threads=${THREADS} \
    --chromosome=chr${CHROM}

log_message "**** Job ends ****"
