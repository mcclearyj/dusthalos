#!/bin/sh
#SBATCH -t 5:59:59
#SBATCH -N 1
#SBATCH --mem=250G
#SBATCH --partition=short
#SBATCH -J dust_cat_maker
#SBATCH -v
#SBATCH -o slurm-loop-cat.out


###
### Define some environmental variables
###

export PATH='.':$PATH:'/work/mccleary_group/Software/texlive-bin/x86_64-linux'
export PYTHONPATH='.':$PYTHONPATH

export CODEDIR='/work/mccleary_group/dusty_halos/dusthalos_emh'
export CONFIGDIR='/work/mccleary_group/dusty_halos/dusthalos_emh/configs'
export CATALOGDIR='/work/mccleary_group/dusty_halos/catalogs'
export SUBCATDIR="${CATALOGDIR}/prep_cat_sdss/"


echo $PATH
echo $PYTHONPATH
echo $CONFIGDIR


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
echo "Skipping prep_cat_runner for galaxy catalog"
echo ""
#python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_sdss_catalog_config.yaml

for value in {2..5} 
do    
    # Make random catalogs                                                                                                            
    #echo "Running make_random_cat_runner.py for foreground #${value}"
    python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_bg_randoms_config_sdss.yaml

    # Prep the random catalogs                                                                                                        
    echo "Running prep_cat_runner for random galaxy catalogs #${value}"
    python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_sdss_rand_catalog_config.yaml
    mv ${SUBCATDIR}/rand_sdss_bg2_JOINED_catalog.fits ${SUBCATDIR}/rand_sdss_bg2_JOINED_catalog${value}.fits

    # Also move the the catalog
    mv ${CATALOGDIR}/sdss_bg_photoz2_randoms.fits ${CATALOGDIR}/sdss_bg_photoz2_randoms${value}.fits
done

# Concatenate catalogs 
echo "Concatenating random galaxy catalogs"
python $CODEDIR/src/concatenate_catalogs.py

echo All done

###
### Record end time 
###                                                                                                                             
echo "Code end time: "
date "+%Y-%m-%d %H:%M:%S"
