#!/bin/bash
#This file is called submit-script.sh
#SBATCH --time=0-3:00:00       # run time in days-hh:mm:ss
#SBATCH --chdir=/work/hageno
#SBATCH --error=logs/job.%J-%a.err
#SBATCH --output=logs/job.%J-%a.out
#SBATCH --mem-per-cpu=4G         # RAM per cpu, in GB
#SBATCH --cpus-per-task=1        #CPU's per taks 

# Load the desired R env with version number
module load foss/2022b R/4.2.2

# Parse input arguments
args=("$@")
# print out to check if all arguments are properly parsed
echo "Start running: ${args[@]:0:8}"

# Run the simulation with Rscript command
Rscript --vanilla ./code/2_scripts/cluster/summary_stats/run_combine_summary_stats.R "M0" "perms_disp_comp_2000_M0"
Rscript --vanilla ./code/2_scripts/cluster/summary_stats/run_combine_summary_stats.R "M1" "perms_disp_comp_2000_M1"
Rscript --vanilla ./code/2_scripts/cluster/summary_stats/run_combine_summary_stats.R "M2" "perms_disp_comp_2000_M2"
Rscript --vanilla ./code/2_scripts/cluster/summary_stats/run_combine_summary_stats.R "M1" "perms_disp_comp_st_2000_M1"