#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=flare_model
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=24gb
#SBATCH --output=logs/flare.%A-%a.log
#SBATCH --array=1-22
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

module purge
module load bcftools/1.10.1
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

mkdir -p "$OUTDIR"
mkdir -p "$TEMPDIR"

log_message "**** Match Variants ****"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
         "${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz" > \
         "${TEMPDIR}/ref.chr${CHROM}.variants.txt"

bcftools view -R "${TEMPDIR}/ref.chr${CHROM}.variants.txt" \
         "${INPUTS}/chr${CHROM}.vcf.gz" | \
    bcftools sort | \
    bcftools norm -d both -Oz -o "${INPUTS}/chr${CHROM}.filtered.vcf.gz"

bcftools view -R "${TEMPDIR}/ref.chr${CHROM}.variants.txt" \
         "${TEMPDIR}/1kGP.chr${CHROM}.filtered.snpsOnly.afr_washington.vcf.gz" | \
    bcftools sort | \
    bcftools norm -d both -Oz -o "${TEMPDIR}/ref.chr${CHROM}.filtered.vcf.gz"

log_message "**** Index VCF files ****"
tabix -f -p vcf "${TEMPDIR}/ref.chr${CHROM}.filtered.vcf.gz"
tabix -f -p vcf "${INPUTS}/chr${CHROM}.filtered.vcf.gz"

log_message "**** FLARE Local Ancestry Analysis ****"

java -Xmx24g -jar $SOFTWARE/flare.jar \
     ref="${TEMPDIR}/ref.chr${CHROM}.filtered.vcf.gz" \
     ref-panel="${TEMPDIR}/samples_id2" nthreads=$THREADS \
     map="${MAPDIR}/plink.chr${CHROM}.GRCh38.reformatted.map" \
     gt="${INPUTS}/chr${CHROM}.filtered.vcf.gz" em=false gen=10 \
     seed=13131313 array=true out="${OUTDIR}/chr${CHROM}"

log_message "**** Job ends ****"
