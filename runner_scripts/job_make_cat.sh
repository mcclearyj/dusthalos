#!/bin/sh
#SBATCH -t 6:00:00
#SBATCH -N 1
#SBATCH --mem=50G
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

conda activate dusthalos


###
### Go!
###

python $CODEDIR/runner_scripts/cat_prep_runner.py --config_path "$CONFIGDIR/prep_hiz_randoms_config.yaml"

mv slurm-dust_cat_maker.out "$dirname"

