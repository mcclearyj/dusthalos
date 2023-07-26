import numpy as np
import healpy as hp
import astropy.units as u
from astropy.io import fits
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.table import Table
import matplotlib.pyplot as plt; plt.ion()
import os

# Define the NSIDE values for the masks
nside_mask_256 = 256
nside_mask_4096 = 4096

# Load the mask files (assuming you have filenames for the masks)
catpath = '/Users/j.mccleary/Research/dusty_halos/catalogs'
bg_mask = 'y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits'
fg_mask = 'WISExSCOSmask.fits'
mask_256 = hp.read_map(os.path.join(catpath, fg_mask))
mask_4096 = hp.read_map(os.path.join(catpath, bg_mask), nest=False, partial=True)

# How many random points? Try to get approximately nrand, if possible.
nrand = int(1e6)
fcover = np.sum(mask_4096 > 0)*1./mask_4096.size
ndraw = np.ceil((nrand/fcover)*1.2).astype(int)

# Generate random pixel indices for the highest resolution (NSIDE=4096)
rng = np.random.default_rng(seed=191477)
random_pixels_4096 = rng.integers(0, hp.nside2npix(nside_mask_4096), 2*ndraw)

# Get the ICRS-frame RA and Dec coordinates of the random pixels for NSIDE=4096
ra_4096, dec_4096 = hp.pix2ang(nside_mask_4096,
                        random_pixels_4096, lonlat=True, nest=False)

# Convert the NSIDE=4096 coordinates to SkyCoord object in ICRS frame
icrs_coords = coord.SkyCoord(ra=ra_4096 * u.deg,
                                dec=dec_4096 * u.deg, frame='icrs')

# Convert the SkyCoord object to Galactic coordinates
galactic_coords = icrs_coords.transform_to('galactic')

# Get the Galactic RA and Dec values
ra_galactic = galactic_coords.l.deg
dec_galactic = galactic_coords.b.deg

# Convert the Galactic coordinates to NSIDE=256 pixels
ipix_256 = hp.ang2pix(nside_mask_256, ra_galactic,
                        dec_galactic, lonlat=True, nest=False)

# Check if the random positions fall within the masked region for both resolutions
ra = []; dec = []
for i in range(2*ndraw):
    if (mask_256[ipix_256[i]]>0) and (mask_4096[random_pixels_4096[i]]>0):
        ra.append(ra_4096[i])
        dec.append(dec_4096[i])

# Make it a SkyCoord object and save it to file
fgr_sky = SkyCoord(ra, dec, frame='icrs', unit=u.deg)
fgr_tab = Table(fgr_sky.to_table())
fgr_tab.write('WISExSCOS_fg_randoms.fits', format='fits', overwrite=True)

# Finally, do a nice overlap plot
fig, ax = plt.subplots(1,1, figsize=(11.5, 8), tight_layout=True, \
                subplot_kw=dict(projection='aitoff'))
ax.grid(True)
ax.set_xlabel('RA', fontsize=14); ax.set_ylabel('Dec', fontsize=14)

# Plot the points - it takes a long time for them all to show up!
ax.plot(fgr_sky.ra.wrap_at('180d').radian, fgr_sky.dec.radian, '.', \
            markersize=0.05,label='Foreground randoms', color='xkcd:neon red')
ax.legend(markerscale=200, loc='upper right', fontsize=14)
fig.tight_layout()
fig.savefig('WISExSCOS_fg_randoms.png')
