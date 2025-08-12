#!/bin/bash
#SBATCH --account=p32505
#SBATCH --partition=gengpu
#$BATCH --gres=gpu:h100:1
#SBATCH --job-name=gen_binaries
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.jr1@northwestern.edu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25gb
#SBATCH --output=log_files/binaries.log
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

module load python/3.10.1
module load gcc/11.2.0
module load cuda/11.6.2-gcc-12.3.0

module list

# Set path variables
VENV_PATH="/projects/p32505/opt/venv/.gpu-dev"

# Activate virtual environment
log_message "**** Activate virtual environment ****"
source "${VENV_PATH}/bin/activate"

if [ $? -ne 0 ]; then
    log_message "Error: Failed to activate virtual environment"
    exit 1
fi

log_message "**** Generate binaries ****"
create-binaries ./

if [ $? -ne 0 ]; then
    log_message "Error: Execution failed"
    exit 1
fi

deactivate
log_message "**** Job ends ****"
