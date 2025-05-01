#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#$BATCH --gres=gpu:h100:1
#SBATCH --job-name=viz_global
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --output=logs/global-ancestry.%j.log
#SBATCH --time=02:00:00

# Function to echo with timestamp
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
module load cuda/12.4.1-gcc-12.3.0
module list

# Set path variables
log_message "**** Loading mamba environment ****"
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
eval "$(mamba shell hook --shell bash)"

ENV_PATH="/projects/p32505/opt/env/AI_env"

mamba activate "$ENV_PATH"
python ../_h/01.global_ancestry.py

if [ $? -ne 0 ]; then
    log_message "Error: mamba or script execution failed"
    exit 1
fi

mamba deactivate
log_message "**** Job ends ****"
