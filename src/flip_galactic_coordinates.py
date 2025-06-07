"""
This is a ChatGPT scipt to flip WISExSCOS coordinates to test 
for residual spatial correlations. It's dumb and specific for 
now, but could be generalized without too much trouble.
"""

from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
from astropy.table import Table

# === Load foreground WISE catalog ===
wise = Table.read(
    '/projects/mccleary_group/dusty_halos/catalogs/prep_cat_sdss/DoubleMasked_wiseScosPhotoz160708_zlt0.15_rCal_gt_17.fits',
    format = 'fits', memmap=True
)
l_fg = wise['l'] 
b_fg = wise['b']

# === Convert to Galactic coordinates ===
coords = SkyCoord(l=l_fg*u.deg, b=b_fg*u.deg, frame='galactic')
l = coords.galactic.l.deg  # Galactic longitude
b = coords.galactic.b.deg  # Galactic latitude

# === Flip in Galactic coordinates ===
# Option 1: Flip longitude (180-degree rotation)
l_flipped = (l + 180.0) % 360.0
#b_flipped = b

# Option 2: Flip latitude (reflect about Galactic equator)
#l_flipped = l
b_flipped = -b

# === Convert back to RA/Dec ===
coords_flipped = SkyCoord(l=l_flipped*u.deg, b=b_flipped*u.deg, frame='galactic')

# === Save to file ===
wise.add_columns([l_flipped, b_flipped], names=["l_flipped", "b_flipped"])
wise.meta['comments'].append("added flipped galactic coordinates using flip_galactic_coordinates.py")
wise.meta['comments'].append("use l_flipped with regular b and b_flipped with regular l")
wise.write(
    '/projects/mccleary_group/dusty_halos/catalogs/prep_cat_sdss/flipped_WISExSCOS_zlt0.15_rCal_gt_17.fits'
)
