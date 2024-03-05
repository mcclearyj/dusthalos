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

# First, make the random catalogs
python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_sdss_fg_randoms_config.yaml
python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_sdss_bg_randoms_config.yaml

# Prep the regular catalogs
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_sdss_catalog_config.yaml

# Prep the random catalogs
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_sdss_rand_catalog_config.yaml

# Dust calc
python $CODEDIR/runner_scripts/dust_calc_runner.py -c $CONFIGDIR/dust_calc_config_sdss.yaml
