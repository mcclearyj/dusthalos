import os
from astropy.table import Table, join
import pdb
import numpy as np

catdir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_sdss'
catlist = [
    os.path.join(catdir, 'rand_sdss_bg2_JOINED_catalog.fits'),
    os.path.join(catdir, 'DoubleMasked_sdss_bg_photoz2.fits'),
]

zcut = 0.8

for cat in catlist:
    catalog = Table.read(cat, format='fits', memmap=True)

    try:
        lowz = catalog['redshift'] > zcut 
    except:
        lowz = catalog['redshift_rand'] > zcut

    colorcut = (catalog['i_corr_csfd'] - catalog['z_corr_csfd']) > -0.103
    wg = lowz * colorcut
    catalog[wg].write(cat.replace(".fits", "_colorcut.fits"), overwrite=True)

