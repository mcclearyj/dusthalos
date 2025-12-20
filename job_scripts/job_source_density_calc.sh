#!/bin/sh
#SBATCH -t 5:59:59
#SBATCH -N 1
#SBATCH --mem=50G
#SBATCH --partition=short
#SBATCH -J source_density_calc
#SBATCH -v
#SBATCH -o slurm-source_density_calc_%j.out


###
### A way to convert time in seconds to a neatly printed day/hour/minute/second format
###
function displaytime {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  (( $D > 0 )) && printf '%d days ' $D
  (( $H > 0 )) && printf '%d hours ' $H
  (( $M > 0 )) && printf '%d minutes ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d seconds\n' $S
}

###
### Activate Conda
###
source /home/j.mccleary/miniconda3/etc/profile.d/conda.sh
conda activate dustyhalos

###
###
### Define some environmental variables
###
export CODEDIR='/projects/mccleary_group/dusty_halos/dusthalos/'
export CONFIGDIR='/projects/mccleary_group/dusty_halos/dusthalos/configs'
export PATH='.':$PATH:'/projects/mccleary_group/Software/texlive-bin/x86_64-linux'
export PYTHONPATH='.':$PYTHONPATH

echo "PATH is set to ${PATH}"
echo "PYTHONPATH is set to ${PYTHONPATH}"
echo "CONFIGDIR is set to ${CONFIGDIR}"

dirname="slurm_outfiles"
if [ ! -d "$dirname" ]
then
     echo " Directory $dirname does not exist. Creating now"
     mkdir -p -- "$dirname"
     echo " $dirname created"
 else
     echo " Directory $dirname exists"
 fi

 echo "Proceeding with code..."
 
### 
### Record start time 
### 
echo "Code start time: "
date "+%Y-%m-%d %H:%M:%S"
StartTime=$(date +%s)

###
### Go!
###

echo "Running source density calculation for redMaGiC hi-z"
python $CODEDIR/runner_scripts/source_density_calc_runner.py -c $CONFIGDIR/source_density_calc_hiz.yaml
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"

echo "Running source density calculation for redMaGiC hi-dens"
python $CODEDIR/runner_scripts/source_density_calc_runner.py -c $CONFIGDIR/source_density_calc_hidens.yaml
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"


###
### Record end time and calculate approximate run time
###  
EndTime=$(date +%s)
Diff=$((EndTime - StartTime))
echo "Total elapsed time:"
displaytime $Diff 

### Move slurm outfile
mv slurm-source_density_calc_%j.out "$dirname"



