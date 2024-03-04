###
### For cases where no HEALPix mask is supplied, create one based on
### RA and Dec of catalog (assuming they are representative, ofc!)
###

from src.radec2healpixels import radec2healpixels

bg_catname = '/Users/j.mccleary/Research/dusty_halos/catalogs/sdss_bg_gals_jemcclear.fits'
fg_catname = '/Users/j.mccleary/Research/dusty_halos/catalogs/sdss_fg_gals_jemcclear.fits'

ra_col = 'ra'
dec_col = 'dec'

radec2healpixels(filename=bg_catname, ra_col=ra_col, dec_col=dec_col)
radec2healpixels(filename=fg_catname, ra_col=ra_col, dec_col=dec_col)
