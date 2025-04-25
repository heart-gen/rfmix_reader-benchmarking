#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=flare_model
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24gb
#SBATCH --output=logs/flare.%A-%a.log
#SBATCH --array=1-22
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
module load htslib/1.16
module list

## Job commands here
TEMPDIR="temp"
INPUTS="vcf-files"
OUTDIR="flare-out"
CHROM=${SLURM_ARRAY_TASK_ID}
THREADS=${SLURM_CPUS_PER_TASK}
SOFTWARE="/projects/p32505/opt/bin"
MAPDIR="/projects/b1213/resources/1kGP/genetic_maps"

mkdir -p $OUTDIR

log_message "**** Index reference VCF files ****"
tabix -f -p vcf ${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz

log_message "**** FLARE Local Ancestry Analysis ****"

java -Xmx24g -jar $SOFTWARE/flare.jar \
     ref="${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz" \
     ref-panel="${TEMPDIR}/samples_id2" nthreads=$THREADS \
     map="${MAPDIR}/plink.chr${CHROM}.GRCh38.reformatted.map" \
     gt="${INPUTS}/chr${CHROM}.vcf.gz" em=false gen=10 \
     seed=13131313 array=true out="${OUTDIR}/chr${CHROM}"

log_message "**** Job ends ****"
