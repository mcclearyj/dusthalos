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

export CODEDIR='/work/mccleary_group/dusty_halos/dusthalos'
export CONFIGDIR='/work/mccleary_group/dusty_halos/dusthalos/configs'
export PATH=$PATH:'/work/mccleary_group/Software/texlive-bin/x86_64-linux'

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
### Go!
###

#python $CODEDIR/runner_scripts/prep_cat_runner.py --config_file "$CONFIGDIR/prep_gaia_catalog_config.yaml"
python $CODEDIR/junk.py
#mv slurm-dust_cat_maker.out "$dirname"
