###
### For cases where no HEALPix mask is supplied, create one based on
### RA and Dec of catalog (assuming these are representative of whole survey)
###

from src.radec2healpixels import radec2healpixels

### Define parameters...
catname = '/work/mccleary_group/dusty_halos/catalogs/sdss_bg_photoz2.fits'
ra_col = 'ra'
dec_col = 'dec'
Nside = 1024

### Run!
radec2healpixels(
    filename=catname, Nside=Nside, ra_col=ra_col, dec_col=dec_col
)
