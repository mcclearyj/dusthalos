#!/bin/sh
#SBATCH -t 12:59:59
#SBATCH -N 1
#SBATCH --mem=150G
#SBATCH --partition=short
#SBATCH -J dust_calc
#SBATCH -v
#SBATCH -o slurm-dust_calc.out


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
### Go!
###

python $CODEDIR/runner_scripts/dust_calc_runner.py -c $CONFIGDIR/dust_calc_systematics.yaml
python $CODEDIR/runner_scripts/dust_calc_runner.py -c $CONFIGDIR/dust_calc_hidens.yaml


### Move slurm outfile
mv slurm-dust_calc.out "$dirname"



