###
### For cases where no HEALPix mask is supplied, create one based on
### RA and Dec of catalog (assuming these are representative of whole survey)
###

from src.radec2healpixels import radec2healpixels

### Define parameters...
catname = '/n23data1/mccleary/catalogs/COSMOSWeb_mastercatalog_v1_clean.fits'
ra_col = 'ra'
dec_col = 'dec'
Nside = 4096

### Run!
radec2healpixels(
    filename=catname, Nside=Nside, ra_col=ra_col, dec_col=dec_col
)
