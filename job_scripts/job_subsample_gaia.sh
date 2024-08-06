#!/bin/sh
#SBATCH -t 59:59
#SBATCH -N 1
#SBATCH --mem=300G
#SBATCH --partition=express
#SBATCH -J gaia_subsample
#SBATCH -v
#SBATCH -o gaia_subsample.out


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


###
### Activate Conda
###

source /work/mccleary_group/miniconda3/etc/profile.d/conda.sh
conda activate dustyhalos

###
### Go!
###

echo "Subsampling Gaia catalog"
python $CODEDIR/runner_scripts/subsample_gaia_cat.py 
