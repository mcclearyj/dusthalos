#!/bin/sh
#SBATCH -t 2:59:59
#SBATCH --mem=150G
#SBATCH --partition=short
#SBATCH -J prep_cat
#SBATCH -v
#SBATCH -o slurm_prep_cat_%j.out

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

# Prep the regular catalog
echo ""
echo "Running prep_cat_runner for galaxy catalog"
echo ""
#python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_gaia_catalog_config.yaml
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_fg_hidens_catalog.yaml

echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"


# Prep the random catalog
echo ""
echo "Running prep_cat_runner for random catalog"
echo ""
#python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_gaia_random_config.yaml
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_fg_hidens_randoms.yaml

echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"


###
### Record end time and total elapsed time
###                                                                                                                           
echo "\n\nCode end time: "
date "+%Y-%m-%d %H:%M:%S"

EndTime=$(date +%s)
Diff=$((EndTime - StartTime))
echo "\n\nTotal elapsed time:"
displaytime $Diff 

###
### Move output file to slurm output directory 
###
mv slurm-prep_cat.out "$dirname"
