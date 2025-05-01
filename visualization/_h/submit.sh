#!/bin/bash

# Set the log file paths
LOG_FILE="gpu_write.log"
ERROR_LOG="gpu_write.error"

# Print some information about the job
echo "Start time: $(date)" > $LOG_FILE

# Activate the virtual environment
source /home/kynon/.venvs/ai_env/bin/activate

# Check if the environment was activated successfully
if [ $? -ne 0 ]; then
    echo "Failed to activate virtual environment. Exiting." >> $ERROR_LOG
    exit 1
fi

echo "Virtual environment 'scrna' activated successfully." >> $LOG_FILE

# Run script
python ../_h/01.write-testing.py >> $LOG_FILE 2>&1

# Check if the script ran successfully
if [ $? -ne 0 ]; then
    echo "Python script failed. Check the error logs." >> $ERROR_LOG
    exit 1
fi

# Deactivate the virtual environment
deactivate

echo "Job finished at: $(date)" >> $LOG_FILE
