#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=flare_model
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16gb
#SBATCH --output=logs/flare.%j.log
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
module load bcftools/1.10.1
module list

## Job commands here
CHROM=1
TEMPDIR="temp"
VCFDIR="vcf-files"
OUTDIR="flare-out"
THREADS=${SLURM_CPUS_PER_TASK}
SOFTWARE="/projects/p32505/opt/bin"
MAP_DIR="/projects/b1213/resources/1kGP/genetic_maps"
REF="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/local-ancestry-ref"

mkdir -p "$OUTDIR"
mkdir -p "$TEMPDIR"

log_message "**** Match Variants ****"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
         "${REF}/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz" > \
         "${TEMPDIR}/ref.chr${CHROM}.variants.txt"

bcftools view -R "${TEMPDIR}/ref.chr${CHROM}.variants.txt" \
         "${VCFDIR}/chr${CHROM}.vcf.gz" | \
    bcftools sort | \
    bcftools norm -d both -Oz -o "${VCFDIR}/chr${CHROM}.filtered.vcf.gz"

bcftools view -R "${TEMPDIR}/ref.chr${CHROM}.variants.txt" \
         "${REF}/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz" | \
    bcftools sort | \
    bcftools norm -d both -Oz -o "${TEMPDIR}/ref.chr${CHROM}.filtered.vcf.gz"

log_message "**** Index VCF files ****"
tabix -f -p vcf "${TEMPDIR}/ref.chr${CHROM}.filtered.vcf.gz"
tabix -f -p vcf "${VCFDIR}/chr${CHROM}.filtered.vcf.gz"

log_message "**** FLARE Local Ancestry Analysis ****"

java -Xmx16g -jar $SOFTWARE/flare.jar \
     ref="${TEMPDIR}/ref.chr${CHROM}.filtered.vcf.gz" \
     ref-panel="$REF/samples_id2" \
     map="${MAPDIR}/plink.chr${CHROM}.GRCh38.reformatted.map" \
     gt="${VCFDIR}/chr${CHROM}.filtered.vcf.gz" nthreads=$THREADS \
     seed=13131313 array=true out="${OUTDIR}/chr${CHROM}"

log_message "**** Job ends ****"
