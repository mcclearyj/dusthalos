#!/bin/sh
#SBATCH -t 5:59:59
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH --partition=short
#SBATCH -J dust_cat_maker
#SBATCH -v
#SBATCH -o slurm-loop-cat.out


###
### Define some environmental variables
###

export CODEDIR='/work/mccleary_group/dusty_halos/dusthalos_emh'
export CONFIGDIR='/work/mccleary_group/dusty_halos/dusthalos_emh/configs'
export CATALOGDIR='/work/mccleary_group/dusty_halos/catalogs'

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

for value in {1..5} 
do    
    # Make random catalogs                                                                                                            
    #echo "Running make_random_cat_runner.py for foreground #${value}"
    python $CODEDIR/runner_scripts/make_random_cat_runner.py -c $CONFIGDIR/make_sdss_bg_randoms_config.yaml

    # Prep the random catalogs                                                                                                        
    echo "Running prep_cat_runner for random galaxy catalogs #${value}"
    python $CODEDIR/runner_scripts/prep_cat_runner.py -c $CONFIGDIR/prep_sdss_rand_catalog_config.yaml
    mv ${CATALOGDIR}/rand_sdss_bg2_JOINED_catalog.fits ${CATALOGDIR}/rand_sdss_bg_JOINED_catalog${value}.fits

    # Also move the the catalog
    mv ${CATALOGDIR}/sdss_bg_photoz_randoms.fits ${CATALOGDIR}/sdss_bg_photoz_randoms${value}.fits
done

# Concatenate catalogs 
echo "Concatenating random galaxy catalogs"
python $CODEDIR/src/concatenate_catalogs.py

echo All done
