#!/bin/sh
#SBATCH -t 23:59:59
#SBATCH -N 1
#SBATCH --mem=60G
#SBATCH --partition=short
#SBATCH -J dust_cat_maker
#SBATCH -v
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jmac.ftw@gmail.com
#SBATCH -o slurm-dust_cat_maker.out


###
### Define some environmental variables
###

export CODEDIR='/Users/j.mccleary/Research/dusty_halos/dusthalos_emh'
export CONFIGDIR='/Users/j.mccleary/Research/dusty_halos/dusthalos_emh/configs'
#export PATH=$PATH:'/Users/j.mccleary/Software/texlive-bin/x86_64-linux'

echo $PATH
echo $PYTHONPATH
echo $CONFIGDIR


###
### Activate Conda
###
source /Users/j.mccleary/Software/miniconda3/etc/profile.d/conda.sh
conda activate dustyhalos_emh

###
### Go!
###

# Prep the regular catalogs
echo "Running prep_cat_runner for regular galaxy catalogs"
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hidens_catalog_config.yaml

# Make the random catalogs
echo "Running make_random_cat_runner.py for foreground"
python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_fg_randoms_config.yaml

# Prep the random catalogs
echo "Running prep_cat_runner for random galaxy catalogs"
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hidens_randoms_config.yaml

# Dust calc
echo "Running dust calc runner"
python $CODEDIR/runner_scripts/dust_calc_runner.py -c $CONFIGDIR/dust_calc_config.yaml
