#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=normal
#SBATCH --job-name=flare_model
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16gb
#SBATCH --output=logs/flare.%j.log
#SBATCH --time=24:00:00

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

## Job commands here
CHROM=1
THREADS=${SLURM_CPUS_PER_TASK}
OUTPUT_PREFIX="simu_3pop"
SOFTWARE="/projects/p32505/opt/bin"
MAP_DIR="/projects/b1213/resources/1kGP/genetic_maps"

log_message "**** FLARE Local Ancestry Analysis ****"

java -Xmx16g -jar $SOFTWARE/flare.jar \
     ref="./temp/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz" \
     ref-panel="./temp/samples_id2" \
     map="${MAP_DIR}/plink.chr${CHROM}.GRCh38.map" \
     gt="chr${CHROM}.vcf.gz" nthreads=$THREADS \
     seed=13131313 array=true out="${OUTPUT_PREFIX}"

log_message "**** Job ends ****"
