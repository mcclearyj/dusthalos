###
### For cases where no HEALPix mask is supplied, create one based on
### RA and Dec of catalog (assuming they are representative, ofc!)
###

from src.radec2healpixels import radec2healpixels

catname = '/work/mccleary_group/dusty_halos/catalogs/gaia_stars_southern_full.fits'
ra_col = 'ra'
dec_col = 'dec'

radec2healpixels(filename=catname, ra_col=ra_col, dec_col=dec_col)

