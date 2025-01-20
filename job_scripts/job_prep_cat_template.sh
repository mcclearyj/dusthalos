#!/bin/sh
#SBATCH -t 1:59:59
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH --partition=short
#SBATCH -J prep_cat
#SBATCH -v
#SBATCH -o slurm-prep_cat.out


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

source /work/mccleary_group/miniconda3/etc/profile.d/conda.sh
conda activate dustyhalos

###
### Record start time 
###                                                                                                                             
echo "Code start time: "
date "+%Y-%m-%d %H:%M:%S"

###
### Go!
###

# Prep the regular catalog
echo ""
echo "Running prep_cat_runner for galaxy catalog"
echo ""
#python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_gaia_catalog_config.yaml
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_fg_hidens_catalog.yaml

# Prep the random catalog
echo ""
echo "Running prep_cat_runner for random catalog"
echo ""
#python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_gaia_random_config.yaml
python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_fg_hidens_randoms.yaml

###
### Record end time  
###                                                                                                                           
echo "Code end time: "
date "+%Y-%m-%d %H:%M:%S"

###
### Move output file to slurm output directory 
###
mv slurm-prep_cat.out "$dirname"
