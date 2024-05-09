from astropy.table import Table
import numpy as np
import os 

catdir = '/work/mccleary_group/dusty_halos/catalogs/'
catname = 'gaia_stars_southern_full.fits' 

# Read in the catalog -- it's large!
stars = Table.read(os.path.join(catdir, catname), format='fits', memmap=True)

# Instantiate an RNG to select indices
rng = np.random.default_rng()

# OK here we go. 
# Make 3 sub-samples of 20M stars from Gaia catalog for whatever needs arise.
size = 2e7

for i in range(1, 5):
    # Define name
    outname = os.path.join(catdir, f'gaia_stars_southern_subsample{i}.fits')
    indices = rng.integers(low=0, high=len(stars), size=int(size))
    
    # Add a dummy redshift column, too. 
    zcol = np.zeros(int(size))
    stars[indices].add_column(zcol, name='z')
    
    # Write to file
    stars[indices].write(outname, format='fits', overwrite=True)
    
