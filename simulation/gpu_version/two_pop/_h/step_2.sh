#!/bin/bash
#SBATCH --partition=gpu,caracol
#SBATCH --gpus=1
#SBATCH --job-name=rfmix_gpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=125gb
#SBATCH --output=gpu_usage.no_binaries.log
#SBATCH --time=03:00:00

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

## Run job

python ../_h/01.rfmix_reader.py

echo "**** Job ends ****"
date
