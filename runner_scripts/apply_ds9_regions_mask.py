from regions import Regions
import numpy as np
import healpy as hp
from astropy.io import fits
from astropy.table import Table, vstack, hstack
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import pdb
import os

# For convenience, define file paths here
#catalog_path = '/Users/j.mccleary/Research/jwst_cosmos/catalogs/'
#image_path = '/Users/j.mccleary/Research/jwst_cosmos/real_data/v0.8/'
#image = 'mosaic_nircam_f150w_COSMOS-Web_30mas_A1_v0_8_i2d.fits.gz'
catalog_path = '/n23data1/mccleary/catalogs/'
image_path = '/n23data2/cosmosweb-public/DR1/data/nircam'
image = 'mosaic_nircam_f277w_COSMOS-Web_30mas_A7_v1.0_i2d.fits.gz'

# First load regions file
regions = Regions.read(os.path.join(catalog_path, 'MASK_HSC-stars_inside_CWfootprint.reg'))

# Load image WCS
hdr = fits.getheader(os.path.join(image_path, image), ext=1)
w = WCS(hdr)

# Load random catalog and get coordinates
random_cat = Table.read(os.path.join(catalog_path, 'clean_mastercat_v1_randoms-try1.fits'))
random_coords = SkyCoord(ra=random_cat['ra'], dec=random_cat['dec'], unit='deg', frame='icrs')
random_mask = np.ones(len(random_cat), dtype=bool)

# Loop through regions and mask random catalog
for region in regions:
    mask = region.contains(random_coords, w) # Get indices of points inside region
    if True in mask:
        random_mask = random_mask & ~mask
        print(f'Masked {np.sum(mask)} random points for region: {region}')

# Apply mask to random catalog
masked_randoms = random_cat[random_mask]
masked_randoms.write(
    os.path.join(catalog_path, 'clean_mastercat_v1_randoms-try1_masked.fits'),
    overwrite=True
)