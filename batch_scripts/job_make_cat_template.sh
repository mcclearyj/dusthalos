#!/bin/sh
#SBATCH -t 59:59
#SBATCH -N 1
#SBATCH --mem=120G
#SBATCH --partition=express
#SBATCH -J overlap_plotter
#SBATCH -v
#SBATCH -o slurm-cat_maker.out


###
### Define some environmental variables
###


export CODEDIR='/work/mccleary_group/dusty_halos/dusthalos_emh'
export CONFIGDIR='/work/mccleary_group/dusty_halos/dusthalos_emh/configs'
export PATH='.':$PATH:'/work/mccleary_group/Software/texlive-bin/x86_64-linux'
export PYTHONPATH='.':$PYTHONPATH

echo $PATH
echo $PYTHONPATH
echo $CONFIGDIR


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
### Activate Conda
###

#source /Users/j.mccleary/Software/miniconda3/etc/profile.d/conda.sh
source /work/mccleary_group/miniconda3/etc/profile.d/conda.sh
conda activate dustyhalos

###
### Go!
###

echo "Running prep_cat_runner for real galaxy catalogs"
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hiz_catalog_config.yaml

echo "Plotting overlap plots" 
python $CODEDIR/runner_scripts/make_overlap_plots.py

mv slurm-dust_cat_maker.out "$dirname"
