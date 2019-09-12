#!/bin/bash
#PBS -q mediumq
#PBS -l select=1:ncpus=16
#PBS -j oe
#PBS -o dustcorr.out

# NOTES
# '#PBS' directives must immediately follow your shell initialization line '#!/bin/<shell>'
# '#PBS' directives must be consecutively listed without any empty lines in-between directives
# Please reference the PBS Pro User Guide for #PBS directive options or type ‘man qsub’ on
#  the login node.
# To determine the appropriate queue (#PBS -q) and walltime (#PBS -l walltime) to use,
#  run (qmgr -c 'print server') on the login node.

# In this example, the job is submitted to the shortq, requests 1 node with 8 cores,
#  and a max wall clock time of 1 hour
# By default, stdout and stderr are written to separate files with the job ID in the
#  file name. You can optionally join the stdout and stderr (#PBS -j oe) and name the
#  output file (#PBS -o <name_of_file>) 

# This is an example job script for the JPL Aurora/Halo cluster

# Set the output directory
# By default, PBS copies stdout and stderr files back to $PBS_O_WORKDIR
# When 'qsub' is run, PBS sets $PBS_O_WORKDIR to the directory where qsub is run.
# Change this environment variable if desired
#

export PBS_O_WORKDIR=/home/jemcclea/data2/des_dust/dusthalos/src
# Set your executable directory (optional)
#
export RUN_DIR=/home/jemcclea/data2/des_dust/dusthalos/src

# Load software modules
# Available modules can be found with the command 'module avail'
# You must first load the appropriate init script to load modules environment variable
# Module init scripts are located at /usr/share/Modules/init/ on Aurora and Halo
#
source /usr/share/Modules/init/bash
source /home/jemcclea/.bashrc
#module load compilers/intel-13.1.3

#Run your application
echo $RUN_DIR
cd $RUN_DIR

python dust_measurement.py
