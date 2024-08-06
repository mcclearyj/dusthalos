import os
from astropy.table import Table, join
import pdb
import numpy as np

catdir = '/work/mccleary_group/dusty_halos/catalogs/'
catlist = [
    os.path.join(catdir, 'y3a2_gold2.2.1_redmagic_highdens.fits'),
    os.path.join(catdir, 'y3a2_gold2.2.1_redmagic_highdens_randoms.fits'),
    os.path.join(catdir, 'prep_cat_csfd/redmagic_hidens_randoms_y3_GOLD_JOINED_catalog.fits'), 
    os.path.join(catdir,'prep_cat_csfd/redmagic_hidens_y3_GOLD_JOINED_catalog.fits')
]

zcut = 0.45

for cat in catlist:
    catalog = Table.read(cat, format='fits', memmap=True)
    try:
        lowz = catalog['zredmagic'] < zcut 
        hiz = catalog['zredmagic'] > (zcut + 0.05)
    except:
        lowz = catalog['z'] < zcut 
        hiz = catalog['z'] > (zcut + 0.05)    
    catalog[lowz].write(cat.replace(".fits", "_z_lt_0.45.fits"), overwrite=True)
    catalog[hiz].write(cat.replace(".fits", "_z_gt_0.5.fits"), overwrite=True)

