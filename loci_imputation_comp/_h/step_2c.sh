#!/bin/bash
#SBATCH --partition=RM-shared
#SBATCH --job-name=phase_merged_data
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kj.benjamin90@gmail.com
#SBATCH --ntasks-per-node=64
#SBATCH --time=12:00:00
#SBATCH --output=logs/merge_phased.%j.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

log_message "**** Bridges-2 info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURM_NODENAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

module purge
module load anaconda3/2024.10-1
module list

log_message "**** Loading conda environment ****"
conda activate /ocean/projects/bio250020p/shared/opt/env/ai_env

log_message "**** Run analysis ****"
RFMIX_TWO="input/simulations/two_populations/_m/rfmix-out/phased_files"
RFMIX_THREE="input/simulations/three_populations/_m/rfmix-files/phased_files"

merge-phased-zarrs "${RFMIX_TWO}/phased_data.zarr" \
		   "${RFMIX_TWO}/phased_chr1.zarr" "${RFMIX_TWO}/phased_chr2.zarr" \
		   "${RFMIX_TWO}/phased_chr3.zarr" "${RFMIX_TWO}/phased_chr4.zarr" \
		   "${RFMIX_TWO}/phased_chr5.zarr" "${RFMIX_TWO}/phased_chr6.zarr" \
		   "${RFMIX_TWO}/phased_chr7.zarr" "${RFMIX_TWO}/phased_chr8.zarr" \
		   "${RFMIX_TWO}/phased_chr9.zarr" "${RFMIX_TWO}/phased_chr10.zarr" \
		   "${RFMIX_TWO}/phased_chr11.zarr" "${RFMIX_TWO}/phased_chr12.zarr" \
		   "${RFMIX_TWO}/phased_chr13.zarr" "${RFMIX_TWO}/phased_chr14.zarr" \
		   "${RFMIX_TWO}/phased_chr15.zarr" "${RFMIX_TWO}/phased_chr16.zarr" \
		   "${RFMIX_TWO}/phased_chr17.zarr" "${RFMIX_TWO}/phased_chr18.zarr" \
		   "${RFMIX_TWO}/phased_chr19.zarr" "${RFMIX_TWO}/phased_chr20.zarr" \
		   "${RFMIX_TWO}/phased_chr21.zarr" "${RFMIX_TWO}/phased_chr22.zarr"

merge-phased-zarrs "${RFMIX_THREE}/phased_data.zarr" \
		   "${RFMIX_THREE}/phased_chr1.zarr" "${RFMIX_THREE}/phased_chr2.zarr" \
		   "${RFMIX_THREE}/phased_chr3.zarr" "${RFMIX_THREE}/phased_chr4.zarr" \
		   "${RFMIX_THREE}/phased_chr5.zarr" "${RFMIX_THREE}/phased_chr6.zarr" \
		   "${RFMIX_THREE}/phased_chr7.zarr" "${RFMIX_THREE}/phased_chr8.zarr" \
		   "${RFMIX_THREE}/phased_chr9.zarr" "${RFMIX_THREE}/phased_chr10.zarr" \
		   "${RFMIX_THREE}/phased_chr11.zarr" "${RFMIX_THREE}/phased_chr12.zarr" \
		   "${RFMIX_THREE}/phased_chr13.zarr" "${RFMIX_THREE}/phased_chr14.zarr" \
		   "${RFMIX_THREE}/phased_chr15.zarr" "${RFMIX_THREE}/phased_chr16.zarr" \
		   "${RFMIX_THREE}/phased_chr17.zarr" "${RFMIX_THREE}/phased_chr18.zarr" \
		   "${RFMIX_THREE}/phased_chr19.zarr" "${RFMIX_THREE}/phased_chr20.zarr" \
		   "${RFMIX_THREE}/phased_chr21.zarr" "${RFMIX_THREE}/phased_chr22.zarr"

if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs."
    exit 1
fi

conda deactivate
log_message "Job finished at: $(date)"
