###
### Manually de-redden magnitudes using the corrected SFD of XXX ref needed XXX
###

from astropy.io import fits
from dustmaps.csfd import CSFDQuery
from astropy.coordinates import SkyCoord
import astropy.units as u
from dust_extinction.parameter_averages import G23
import numpy as np

def deredden(catname, wavelengths):
    """ Cat should be a string object """

    # Open catalog, grab stuff 
    cat = fits.open(catname, mode='update', save_backup=True)
    ra = cat[1].data['ra']
    dec = cat[1].data['dec']
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

    # First, define bandpass column names in catalog
    band_names = [
        'mof_cm_mag_g', 'mof_cm_mag_r',
        'mof_cm_mag_i', 'mof_cm_mag_z'
    ]

    # DES also has some bonus ones
    delta_mag_y4_names = [
        'delta_mag_y4_g', 'delta_mag_y4_r',
        'delta_mag_y4_i', 'delta_mag_y4_z'
    ]
    delta_mag_chrom_names = [
        'delta_mag_chrom_g', 'delta_mag_chrom_r',
        'delta_mag_chrom_i', 'delta_mag_chrom_z'
    ]

    # Let's go through band by band and add these column names
    # Initialize a dict
    print("Beginning magnitude corrections...")
    corr_band_dict = {}

    for i in range(len(band_names)):

        # Do actual correction
        chromcorr = data[delta_mag_y4_names[i]] + data[delta_mag_chrom_names[i]]
        this_Ax = Ax[i,:]
        corr_band = data[band_names[i]] + chromcorr - this_Ax

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

def main():
    # Central wavelengths (spectrum) to model; should match number of bands!
    wavelengths = [4796.6, 6382.6, 7769.0, 9108.2] * u.AA    

    cat1 = "/work/mccleary_group/dusty_halos/catalogs/18176_desy3_gold_gals.fits"
    deredden(cat1, wavelengths)
    print(f"Finished working on {cat1}")

    # Do one more
    cat2 = "/work/mccleary_group/dusty_halos/catalogs/18175_des_stars.fits"
    deredden(cat2, wavelengths)
    print(f"Finished working on {cat2}")

if __name__ == '__main__':
    main()
    print("All done!")