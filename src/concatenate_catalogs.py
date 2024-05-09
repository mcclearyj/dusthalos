import pdb
import glob
from astropy.table import Table, vstack

catalog_dir = "/work/mccleary_group/dusty_halos/catalogs"

randoms = glob.glob(os.path.join(
    catalog_dir, "sdss_bg_photoz_randoms*.fits"
))
masked_randoms = glob.glob(os.path.join(
    catalog_dir, "prep_cat_output/rand_sdss_bg_JOINED_catalog*.fits"
))

# First, concatenate randoms
holder = []
for rand in randoms:
    holder.append(Table.read(rand, format='fits', memmap=True))

vstack(holder).write(
    os.path.join(catalog_dir, "sdss_bg_photoz_randoms_stacked.fits"), 
    format="fits", overwrite=True
)

# Then, concatenate masked_randoms
holder = []
for rand in masked_randoms:
    holder.append(Table.read(rand, format='fits', memmap=True))

vstack(holder).write(
    os.path.join(catalog_dir, "prep_cat_output/rand_sdss_bg_JOINED_catalog_stacked.fits")
    format="fits", overwrite=True
)


