#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=prep_samples
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --output=logs/samples_prep.log
#SBATCH --time=00:05:00

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
REF="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/local-ancestry-ref"
INPUTS="simulation-files"
OUTPUT="temp"

log_message "**** Prepare samples ****"

mkdir -p ${OUTPUT}

grep -E 'CEU|YRI|PUR' ${INPUTS}/1k_sampleinfo.tsv > ${OUTPUT}/samples_id2
cut -f1 ${OUTPUT}/samples_id2 > ${OUTPUT}/samples_id
cp -v $REF/genetic_map38 ${OUTPUT}/

log_message "**** Job ends ****"
