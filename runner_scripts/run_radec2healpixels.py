###
### For cases where no HEALPix mask is supplied, create one based on
### RA and Dec of catalog (assuming these are representative of whole survey)
###

from src.radec2healpixels import radec2healpixels

### Define parameters...
catname = '/work/mccleary_group/dusty_halos/catalogs/sf_gswlc_galaxies.fits'
ra_col = 'RA'
dec_col = 'Dec'
Nside = 2048

### Run!
radec2healpixels(
    filename=catname, Nside=Nside, ra_col=ra_col, dec_col=dec_col
)

catname = '/work/mccleary_group/dusty_halos/catalogs/quiescent_gswlc_galaxies.fits'

### Run!
radec2healpixels(
    filename=catname, Nside=Nside, ra_col=ra_col, dec_col=dec_col
)
