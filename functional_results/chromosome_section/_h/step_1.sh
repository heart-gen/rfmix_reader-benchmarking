#!/bin/bash
#SBATCH --partition=gpu,caracol
#SBATCH --job-name=visualization
#SBATCH --mail-type=FAIL
#SBATCH --gpus=1
#SBATCH --mail-user=jbenja13@jh.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25gb
#SBATCH --output=summary.log
#SBATCH --time=01:00:00

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

echo "**** Run visualization ****"

cp -v ../_h/01.visualization.ipynb .
quarto render 01.visualization.ipynb --execute
rm 01.visualization.ipynb

echo "**** Job ends ****"
date
