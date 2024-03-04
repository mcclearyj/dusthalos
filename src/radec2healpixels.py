import numpy as np
import healpy as hp
from astropy.table import Table

def radec2healpixels(filename, healpix_name=None,
                     ra_col='ra', dec_col='dec', hdu=1):
    """
    Make a HEALPix mask for a catalog using lists of RA, Dec.
    """
    # Open file
    tab = Table.read(filename, hdu=hdu, memmap=True)

    # Example RA and Dec data (in degrees)
    ra = tab[ra_col]  # Replace with your RA data
    dec = tab[dec_col]  # Replace with your Dec data

    # Convert RA and Dec to radians
    ra_rad = np.radians(ra)
    dec_rad = np.radians(dec)

    # Convert RA, Dec to theta (colatitude) and phi (longitude)
    theta = np.pi/2 - dec_rad  # theta = 90 degrees - Dec
    phi = ra_rad

    # Choose a resolution (Nside)
    Nside = 64  # for example, can be adjusted based on your needs

    # Convert to HEALPix pixel indices for your RA, Dec coordinates
    pixels = hp.ang2pix(Nside, theta, phi)

    # Initialize a HEALPix map with zeros
    healpix_map = np.zeros(hp.nside2npix(Nside), dtype=int)

    # Set the pixels in your list to 1
    for pix in pixels:
        healpix_map[pix] = 1

    if healpix_name == None:
        #healpix_name = "GSWLC-A2_healpix_mask.fits"
        healpix_name = filename.replace('.fits', '_healpix_mask.fits')

    # Save the map to a FITS file
    hp.write_map(healpix_name, healpix_map, overwrite=True)
