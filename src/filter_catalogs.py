import os
from astropy.table import Table, join
import pdb
import numpy as np

catdir = '/work/mccleary_group/dusty_halos/catalogs/'
catlist = [
    os.path.join(catdir, 'sdss_bg_photoz2.fits'),
]

for cat in catlist:
    catalog = Table.read(cat, format='fits', memmap=True)

    #try:
    #    lowz = catalog['redshift'] > zcut 
    #except:
    #    lowz = catalog['redshift_rand'] > zcut

    colorcut = (catalog['i_corr_csfd'] - catalog['z_corr_csfd']) > -0.103
    #wg = lowz * colorcut
    wg = colorcut 
    catalog[wg].write(cat.replace(".fits", "_colorcut.fits"), overwrite=True)

