#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=short
#SBATCH --job-name=flare_prep
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.jr1@northwestern.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=16gb
#SBATCH --output=log_files/flare.%A-%a.log
#SBATCH --array=2-22
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
OUTDIR="flare-out"
CHROM=${SLURM_ARRAY_TASK_ID}
THREADS=${SLURM_CPUS_PER_TASK}
SOFTWARE="/projects/p32505/opt/bin"
TEMP_MAP_DIR="/projects/p32505/users/manuel/rfmix_reader-benchmarking/input/simulations/two_populations/_m/temp"
REF="/projects/b1213/resources/1kGP/GRCh38_phased_vcf/local-ancestry-ref"

log_message "**** FLARE Local Ancestry Analysis ****"

java -Xmx16g -jar $SOFTWARE/flare.jar \
     ref="$REF/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz" \
     ref-panel="$REF/samples_id2" \
     map="$TEMP_MAP_DIR/plink.chr${CHROM}.GRCh38.reformatted.map" \
     gt="chr${CHROM}.vcf.gz" nthreads=$THREADS \
     seed=13131313 em=false model="${OUTDIR}/chr1.model" \
     array=true out="${OUTDIR}/chr${CHROM}"

log_message "**** Job ends ****"
