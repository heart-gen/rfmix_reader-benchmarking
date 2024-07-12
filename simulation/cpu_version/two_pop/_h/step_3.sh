#!/bin/bash
#SBATCH --partition=bluejay,shared
#SBATCH --job-name=pandas_cpu
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=450gb
#SBATCH --output=cpu_usage.pandas.log
#SBATCH --time=12:00:00

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

python ../_h/02.pandas.py

echo "**** Job ends ****"
date
