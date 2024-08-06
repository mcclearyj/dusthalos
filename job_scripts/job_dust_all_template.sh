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
echo ""
echo "Running prep_cat_runner for regular galaxy catalogs"
echo ""
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hiz_catalog_config.yaml

# Make the random catalogs
echo ""
echo "Running make_random_cat_runner.py for foreground"
echo ""
python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_fg_randoms_config.yaml

# Prep the random catalogs and do a comfort plot
echo ""
echo "Running prep_cat_runner for random galaxy catalogs"
echo ""
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_hiz_randoms_config.yaml

echo ""
echo "Make overlap plots"
echo ""
python $CODEDIR/runner_scripts/make_overlap_plots_local.py

# Dust calculation(s)
echo ""
echo "Running dust calc runner: hidens"
echo ""
python $CODEDIR/runner_scripts/dust_calc_runner.py -c $CONFIGDIR/dust_calc_config.yaml


