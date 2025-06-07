###
### Manually de-redden magnitudes using the corrected SFD of XXX ref needed XXX
###

from astropy.io import fits
from dustmaps.csfd import CSFDQuery
from astropy.coordinates import SkyCoord
import astropy.units as u
from dust_extinction.parameter_averages import G23
import numpy as np

def deredden(**params):
    """ Parse arguments, run deredden """

    try:
        catname = params['catalog']
        band_names = params['band_names']
        wavelengths = params['wavelengths']*u.AA
        ra_colname = params['ra_colname']
        dec_colname = params['dec_colname']

        _deredden(catname, band_names, wavelengths, ra_colname, dec_colname)

    except KeyError as e:
        print("deredden_magnitudes missing required keyword, double-check call")
        print(e)

def _deredden(catname, band_names, wavelengths, ra_colname, dec_colname):
    # Open catalog, grab stuff
    cat = fits.open(catname, mode='update', save_backup=True)
    ra = cat[1].data[ra_colname]
    dec = cat[1].data[dec_colname]
    coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
    print("Got catalog, coordinates")

    # Get our E(B-V) values
    csfd = CSFDQuery()
    ebv = csfd(coords)
    print("Queried coordinates")

    # Convert to Av values; second line makes it a column vector for broadcasting!
    Rv = 3.1
    Av_values = ebv * Rv

    # Do dust modeling
    gordon23 = G23(Rv=Rv)
    AxAv = gordon23.evaluate(wavelengths, Rv)

    # Finally, calculate offsets. np.newaxis makes it a column vector
    Ax = Av_values[:, np.newaxis] * AxAv
    Ax = Ax.T # There is a smarter way to do this, but I don't know it.
    print("Obtained Ax for catalog")

    # OK, time to do corrected magnitudes
    # For convenience
    data = cat[1].data

    # Let's go through band by band and add these column names
    # Initialize a dict
    print("Beginning magnitude corrections...")
    corr_band_dict = {}

    for i in range(len(band_names)):

        # Do actual correction
        this_Ax = Ax[i,:]
        corr_band = data[band_names[i]] - this_Ax

        # Define a key name, extend dict with it
        key_name = f'{band_names[i]}_corr_csfd'
        corr_band_dict[key_name] = corr_band

    # Now make this an HDU, write to file
    print(f"Adding new magnitude columns to catalog {catname}")
    columns = data.columns
    for key, value in corr_band_dict.items():
        try:
            col = fits.Column(name=key,  format='D', array=value)
            columns.add_col(col)
        except ValueError:
            pass

    # Those changes all get saved, as far as I can tell!
    print(f"Saving output")
    cat.flush()
