#!/bin/bash
#SBATCH --partition=gpu,caracol
#SBATCH --gpus=1
#SBATCH --job-name=cudf_gpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb
#SBATCH --output=gpu_usage.cudf.log
#SBATCH --time=01:30:00

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

python ../_h/03.cudf.py

echo "**** Job ends ****"
date
