#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --account=def-fmaps
#SBATCH --job-name=mitgcmuv_all
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=maxime.benoit-gagne@takuvik.ulaval.ca
#SBATCH --array=1

# run in
# /home/benoitga/projects/def-fmaps/benoitga/gud_groups/gud_1d_35+16/input_noradtrans

echo "Starting task ${SLURM_ARRAY_TASK_ID}"

# set constants
GUD_FILES="gud_files_exp0.txt"
RUN_DIRS="run_dirs_exp0.txt"

# test current directory
if [[ $PWD != */input_noradtrans ]]
then
    echo 'ERROR: This file should be run in input_noradtrans'
    exit 1
fi

# run mitgcmuv
run_dir=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${RUN_DIRS})
gud_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${GUD_FILES})
./mitgcmuv_one_exp0.sh "${run_dir}" "${gud_file}"
