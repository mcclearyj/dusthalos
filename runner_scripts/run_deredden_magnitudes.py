###
### Apply extinction correction to generic catalog magnitudes
###

from src.deredden_magnitudes import deredden

### Define input parameters ###
# Catalog to deredden
catname = "/Users/j.mccleary/Research/dusty_halos/catalogs/sdss_fg_gals.fits"
# First, define bandpass column names in catalog
band_names = ['u', 'g', 'r', 'i', 'z']
# Central wavelengths (spectrum) to model; should match number of bands!
wavelengths = [3586.8, 4716.7, 6165.1, 7475.9, 8922.9]
# RA column name
ra_colname = 'ra'
# Dec column name
dec_colname = 'dec'

deredden_args = {
'catname': catname,
'band_names': band_names,
'wavelengths': wavelengths,
'ra_colname': ra_colname,
'dec_colname': dec_colname
}

### Run! ###
deredden(deredden_args)

print(f"Finished working on {catalog}!")
