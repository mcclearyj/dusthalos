#!/bin/sh
#SBATCH -t 6:59:59
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=250G
#SBATCH --partition=short
#SBATCH -J dust_calc_all
#SBATCH -v
#SBATCH -o slurm-dust_calc_all.out


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

# Prep the regular catalogs
echo ""
echo "Running prep_cat_runner for regular galaxy catalogs"
echo ""
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hiz_catalog_config.yaml
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"

# Make the random catalogs
echo ""
echo "Running make_random_cat_runner.py for foreground"
echo ""
python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_fg_randoms_config.yaml
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"

# Prep the random catalogs and do a comfort plot
echo ""
echo "Running prep_cat_runner for random galaxy catalogs"
echo ""
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hiz_randoms_config.yaml
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"

echo ""
echo "Make overlap plots"
echo ""
python $CODEDIR/runner_scripts/make_overlap_plots_local.py
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"

# Dust calculation(s)
echo ""
echo "Running dust calc runner: hidens"
echo ""
python $CODEDIR/runner_scripts/dust_calc_runner.py -c $CONFIGDIR/dust_calc_config.yaml
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"


###
### Record end time and calculate approximate run time
###  
echo "Task end time: "
date "+%Y-%m-%d %H:%M:%S"

EndTime=$(date +%s)
Diff=$((EndTime - StartTime))
echo "Total elapsed time:"
displaytime $Diff 
