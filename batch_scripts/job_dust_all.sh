#!/bin/sh
#SBATCH -t 12:59:59
#SBATCH -N 1
#SBATCH --mem=150G
#SBATCH --partition=short
#SBATCH -J dust_cat_maker1
#SBATCH -v
#SBATCH -o slurm-dust_calc1.out


###
### Define some environmental variables
###

export CODEDIR='/work/mccleary_group/dusty_halos/dusthalos_emh/'
export CONFIGDIR='/work/mccleary_group/dusty_halos/dusthalos_emh/configs'
export PATH='.':$PATH:'/work/mccleary_group/Software/texlive-bin/x86_64-linux'
export PYTHONPATH='.':$PYTHONPATH

echo $PATH
echo $PYTHONPATH
echo $CONFIGDIR


###
### Activate Conda
###
source /work/mccleary_group/miniconda3/etc/profile.d/conda.sh
conda activate dustyhalos

###
### Go!
###

# Prep the regular catalogs
echo "Running prep_cat_runner for regular galaxy catalogs"
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hiz_catalog_config.yaml

# Make the random catalogs
# echo "Running make_random_cat_runner.py for foreground"
python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_fg_randoms_config.yaml

# Prep the random catalogs
echo "Running prep_cat_runner for random galaxy catalogs"
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hiz_randoms_config.yaml

# Dust calc
#echo "Running dust calc runner"
#python $CODEDIR/runner_scripts/dust_calc_runner.py -c $CONFIGDIR/dust_calc_config.yaml
