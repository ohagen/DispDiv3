#!/bin/bash
#This file is called submit-script.sh
#SBATCH --time=0-24:00:00       # run time in days-hh:mm:ss
#SBATCH --chdir=/work/hageno
#SBATCH --error=logs/job.%J-%a.err
#SBATCH --output=logs/job.%J-%a.out
#SBATCH --mem-per-cpu=4G         # RAM per cpu, in GB
#SBATCH --cpus-per-task=1        #CPU's per taks 

# Load the desired R env with version number
module load foss/2022b R/4.2.2 Zip/3.0
# trick R into using the zip from the module system
export R_ZIPCMD=$(which zip)

# Parse input arguments
args=("$@")
# print out to check if all arguments are properly parsed
echo "Start running: ${args[@]:0:8}"

# Run the simulation with Rscript command
Rscript --vanilla run_combnine_summary_stats "M2" $SLURM_ARRAY_TASK_ID "perms_disp_comp_2000_M2"