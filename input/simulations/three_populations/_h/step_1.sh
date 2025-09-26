#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --job-name=simu_3pops
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --array=1-22
#SBATCH --cpus-per-task=1
#SBATCH --mem=150gb
#SBATCH --output=logs/simgeno.%A_%a.log
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
log_message "**** Loading modules ****"

module purge
module load htslib/1.16
module load bcftools/1.10.1
module list

## Edit with your job command
log_message "**** Loading mamba environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate /projects/p32505/opt/envs/ml

echo "**** Run simulation ****"
ONE_K="/projects/b1213/resources/1kGP/"
CHROM=${SLURM_ARRAY_TASK_ID}
INPUTS="simulation-files"
VCFDIR="vcf-files"
TEMPDIR="temp"

mkdir -p "$VCFDIR"
mkdir -p "$TEMPDIR"

# Symlink map file
ln -s "${ONE_K}/genetic_maps/plink.chr${CHROM}.GRCh38.map" \
   "$TEMPDIR/chr${CHROM}.map"

# Pre-filter the reference VCF to biallelic variants only
FILTERED_VCF="${TEMPDIR}/chr${CHROM}.biallelic.vcf.gz"
bcftools view -m2 -M2 -v snps,indels \
    ${ONE_K}/GRCh38_phased_vcf/raw/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz \
    -Oz -o ${FILTERED_VCF}
bcftools index -f ${FILTERED_VCF}

# Run haptools on biallelic-only VCF
haptools simgenotype \
         --model ${INPUTS}/AFR_washington.dat \
         --mapdir ${TEMPDIR}/ \
         --chroms ${CHROM} \
         --seed 20240126 \
         --ref_vcf ${FILTERED_VCF} \
         --sample_info ${INPUTS}/1k_sampleinfo.tsv \
         --out ${VCFDIR}/chr${CHROM}.vcf.gz

# Index output
tabix -p vcf ${VCFDIR}/chr${CHROM}.vcf.gz

# Move breakpoint file
mv ${VCFDIR}/chr${CHROM}.bp ${INPUTS}/

conda deactivate
log_message "**** Job ends ****"
