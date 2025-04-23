#!/bin/bash
#SBATCH --job-name=rfmix_aaOnly
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --partition=shared,bluejay
#SBATCH --nodes=1
#SBATCH --array=23
#SBATCH --cpus-per-task=1
#SBATCH --mem=25gb
#SBATCH --output=logs/summary.%A_%a.log
#SBATCH --time=72:00:00

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## List current modules for reproducibility

module list

## Edit with your job command

echo "**** Run RFMix ****"
REF="../../inputs/vcf_ref/_m"
INPUT="../../inputs/genotypes/_m/"

if [[ "${SLURM_ARRAY_TASK_ID}" -eq 23 ]]; then
    CHROM="X"
    echo "Chromosome: ${CHROM}"
    rfmix \
	-f $INPUT/TOPMed_LIBD_AA.chr${CHROM}.vcf.gz \
        -r $REF/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz \
	-m $REF/samples_id2 \
	-g $REF/genetic_map38.chrX \
	-o chr$CHROM \
	--chromosome=${CHROM}
else
    CHROM=${SLURM_ARRAY_TASK_ID}
    rfmix \
	-f $INPUT/TOPMed_LIBD_AA.chr${CHROM}.vcf.gz \
	-r $REF/1kGP_high_coverage_Illumina.chr${CHROM}.filtered.SNV_INDEL_SV_phased_panel.snpsOnly.eur.afr.vcf.gz \
	-m $REF/samples_id2 \
	-g $REF/genetic_map38 \
	-o chr$CHROM \
	--chromosome=chr$CHROM
fi

echo "**** Job ends ****"
date -Is
