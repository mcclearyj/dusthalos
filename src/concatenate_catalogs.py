import pdb
import glob
from astropy.table import Table, vstack
import os

catalog_dir = "/work/mccleary_group/dusty_halos/catalogs"

masked_randoms = glob.glob("/scratch/j.mccleary/rand_sdss_bg2_JOINED_catalog[1-3].fits")

# Concatenate masked_randoms
holder = []
for rand in masked_randoms:
    holder.append(Table.read(rand, format='fits', memmap=True))

vstack(holder).write(
    os.path.join(catalog_dir, "prep_cat_sdss/rand_sdss_bg2_JOINED_catalog_stacked.fits"),
    format="fits", overwrite=True
)


