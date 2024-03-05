from astropy.io import fits
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.table import Table
import matplotlib.pyplot as plt; plt.ion()
import numpy as np

# Load in the WISExSCOS catalog with only SCOS mask applied
scos_fg = Table.read('WISExSCOS_fg_SCOSmask_cat.fits')

# Load in the WISExSCOS catalog with redMaGiC mask applied too
redmagic_mask = Table.read('WISExSCOS_SCOSmask_rmzMask_cat.fits')

# Create SkyCoord object to hold RA, Dec of catalogs
gal = SkyCoord(scos_fg['l'], scos_fg['b'], frame='galactic', unit=u.deg)
eq = gal.icrs
sky2 = SkyCoord(redmagic_mask['ra'], redmagic_mask['dec'], frame='icrs', unit=u.deg)

# Create a plot instance (can also use axes class)
# Note: aitoff projection apparently avoids mollweide's extreme edge distortions

fig = plt.subplot(111, projection='aitoff'); plt.grid(True)
plt.xlabel('RA', fontsize=16); plt.ylabel('Dec', fontsize=16)

# Plot the points - it takes a long time for them all to show up!
plt.plot(eq.ra.wrap_at('180d').radian, eq.dec.radian, '.', markersize=0.02, label='SCOSmask only', color='xkcd:bluegrey')
plt.plot(sky2.ra.wrap_at('180d').radian, sky2.dec.radian, '.', markersize=0.02, label='SCOS + redMaGiC masks', color='xkcd:neon red')
plt.legend(markerscale=400, loc='upper right', fontsize=14)
plt.tight_layout()

plt.savefig('WISExSCOS_fg+redmagic_mask.pdf')
plt.savefig('WISExSCOS_fg+redmagic_mask.png')
