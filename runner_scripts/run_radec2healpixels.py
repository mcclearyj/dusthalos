###
### For cases where no HEALPix mask is supplied, create one based on
### RA and Dec of catalog (assuming these are representative of whole survey)
###
import os
from src.radec2healpixels import radec2healpixels

### Define parameters...
catdir = '/projects/mccleary_group/dusty_halos/catalogs'
catname = 'gswlc_quiescent_galaxies.fits'
filename = os.path.join(catdir,catname)
hdu=1
ra_col = 'RA'
dec_col = 'Dec'
Nside = 1024

### Run!
radec2healpixels(
    filename=filename, hdu=hdu,
    ra_col=ra_col, dec_col=dec_col, Nside=Nside
)