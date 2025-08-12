#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=flare_model
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.jr1@northwestern.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16gb
#SBATCH --output=log_files/flare.%j.log
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
CHROM=1
OUTDIR="flare-out"
THREADS=${SLURM_CPUS_PER_TASK}
SOFTWARE="/projects/p32505/opt/bin"
MAP_DIR="/projects/b1213/resources/1kGP/genetic_maps"
REF="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/local-ancestry-ref"

mkdir -p $OUTDIR
mkdir -p temp

log_message "**** Fix PLINK map file ****"
awk '{if(NR>0) $1="chr"$1; print}' "${MAP_DIR}/plink.chr${CHROM}.GRCh38.map" \
    > ./temp/plink.chr${CHROM}.GRCh38.reformatted.map

log_message "**** FLARE Local Ancestry Analysis ****"

java -Xmx16g -jar $SOFTWARE/flare.jar \
     ref="$REF/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz" \
     ref-panel="$REF/samples_id2" \
     map="./temp/plink.chr${CHROM}.GRCh38.reformatted.map" \
     gt="chr${CHROM}.vcf.gz" nthreads=$THREADS \
     seed=13131313 array=true out="${OUTDIR}/chr${CHROM}"

log_message "**** Job ends ****"
