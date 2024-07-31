###
### Manually de-redden magnitudes using the corrected SFD of XXX ref needed XXX
###

from astropy.io import fits
from dustmaps.csfd import CSFDQuery
from astropy.coordinates import SkyCoord
import astropy.units as u
from dust_extinction.parameter_averages import G23
import numpy as np

# Open catalog
cat = fits.open('/Users/j.mccleary/Research/dusty_halos/catalogs/18170_y3gold.fits', mode='update', save_backup=True)
ra = cat[1].data['ra']
dec = cat[1].data['dec']
coords = SkyCoord(ra, dec, frame='icrs', unit='deg')

# Get our E(B-V) values
csfd = CSFDQuery()
ebv = csfd(coords)

# Convert to Av values; second line makes it a column vector for broadcasting!
Rv = 3.1
Av_values = ebv * Rv

# Central wavelengths (spectrum) to model; should match number of bands!
wavelengths = [4808.49, 6417.65, 7814.58, 9168.85] * u.AA

# Do dust modeling --
#waves = [3556.52, 4702.50, 6175.58, 7489.98, 8946.71] *u.AA
gordon23 = G23(Rv=Rv)
AxAv = gordon23.evaluate(wavelengths, Rv)

# Finally, calculate offsets. np.newaxis makes it a column vector
Ax = Av_values[:, np.newaxis] * AxAv
Ax = Ax.T # There is a smarter way to do this, but I don't know it.

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
columns = data.columns
for key, value in corr_band_dict.items():
    try:
        col = fits.Column(name=key,  format='D', array=value)
        columns.add_col(col)
    except ValueError:
        pass

# Those changes all get saved, as far as I can tell!
cat.flush()
