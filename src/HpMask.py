import sys, os
from astropy.table import Table, vstack, hstack
import healpy as hp
import numpy as np
import pdb
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import time
from .cat_utils import mask_config_checker
##
## TO DO: make a single, inheritable seeded random for making random catalogs
## probably using a Mixin class.
##
## In time, I would like to make this a master catalog processor
##

class HpMask:
    '''
    This should hold the essential mask information: mask filename,
    NSIDE, whether mask is partial, the coordinate system in which the mask
    is defined (VERY IMPORTANT!). It should load the mask into memory. Defining
    the empty or unseen HEALPix pixels could be useful too.
    '''

    def __init__(self, filepath=None, partial=False, coordframe=None):

        # Do some sanity checking
        if filepath is None:
            raise ValueError('Missing required parameter "filepath"')
        if coordframe is None:
            print('Assuming mask is based on equatorial coordinates ("icrs")')
            coordframe = 'icrs'

        # Define attributes
        self.filepath = filepath # relative or absolute path
        self.partial = partial
        self.NSIDE = -1
        self.mask = None
        self.mask_header = {}
        self.seen = []
        self.all_nside_hpix = []
        self.coordframe = coordframe

        # Load mask
        self._load_mask(filepath = self.filepath)

    def _load_mask(self, filepath):
        '''
        Load up the mask file with a way to handle partial masks
        '''

        try:
            mask, h = hp.read_map(self.filepath, h=True, partial=self.partial)
        except FileNotFoundError as fnf:
            print("File not found: ", fnf)
        except ValueError as e:
            print("WARNING: hp.read_map() returned the following error: ", e)
            print("This means hpMask is probably partial; setting partial=True")
            self.partial = True
            mask, h = hp.read_map(self.filepath, h=True, partial=self.partial)

        # Set class attributes
        self.mask = mask
        self.mask_header = dict(h)
        self.NSIDE = self.mask_header['NSIDE']
        self.seen = (mask > 0) & (mask != hp.UNSEEN)
        # Incredibly, mask doesn't contain anything like this.
        self.all_nside_hpix = np.arange(hp.nside2npix(self.NSIDE))

    def apply_mask(self, coords, lonlat=True, vb=True):
        '''
        Apply mask to a set of coordinates. This could probably be made into
        a static method.

        Inputs
            coords: SkyCoords instance to be placed in HEALPIx
            lonlat: Boolean specifying whether or not coordinates are
                    longitude/latitude (degrees). If True, they are. If False,
                    they are healpy longitude and co-latitude (units=radians)
            vb: verbose output [True/False; default True]

        Returns
            ind: indices of successful coordinates
            matched_coords: matched coordinates
        '''

        # Make sure coordinates are in acceptable format
        if type(coords) is not SkyCoord:
            raise TypeError('Supplied "coords" must be an instance of ' + \
                            'astropy.coordinates.SkyCoord')

        # OH! I need to do the silly range(len) thing because after being
        # read in by HEALPy,  the mask is just
        # the mask *value* not the mask HEALPixel. So dumb!

        # Identify good mask pixels
        good_map_hpix = self.all_nside_hpix[self.seen]

        # Probably superfluous, but just in case.
        if lonlat is True:
            ra = coords.ra.deg; dec = coords.dec.deg
        else:
            ra = coords.ra.rad; dec = coords.dec.rad

        # Get HEALPixel for each RA, Dec
        hpInd = hp.ang2pix(self.NSIDE, ra, dec, lonlat=lonlat)

        # Identify coordinates that fall into good (unmasked) HEALPixels
        overlap = np.in1d(hpInd, good_map_hpix)

        # Create object index array and return good indices, good coords
        gal_ind = np.arange(len(coords))

        if vb is True:
            print('\n\n HpMask.apply_mask: Mask applied to input SkyCoords')
            print(f" {len(gal_ind[overlap])}/{len(gal_ind)} objects " +
                        "overlapped with mask HEALPix")
            print(f"  --> fractional overlap/match rate = " +
                        f"{len(gal_ind[overlap])*100.0/len(gal_ind):.1f}%" +
                        " of objects \n")

        if len(gal_ind[overlap]) == 0:
            print('WARNING: No galaxies in input catalog matched to mask')

        return gal_ind[overlap]

    @staticmethod
    def apply_overlapping_masks(coords, mask1, mask2):
        '''
        Return coordinates and indices of objects that lie in the overlap area
        of two HEALPix masks. There should probably be an attribute checker
        of some kind...

        Parameters
            coords: SkyCoords instance to be placed in overlapping HEALPixels
            mask1, mask2: Instances of HpMask for which to calculate overlap

        Returns
            good_gals: array indices that fell within an overlapping HEALPixel
        '''

        # Do some sanity checking
        if type(coords) is not SkyCoord:
            raise TypeError('Supplied "coords" is not an instance of ' + \
                            'astropy.coordinates.SkyCoord')
        #if (type(mask1) is not HpMask):
        #    raise TypeError('Supplied "mask1" is not an instance of hpMask')
        #if (type(mask2) is not HpMask):
        #    raise TypeError('Supplied "mask2" is not an instance of hpMask')

        nside1 = mask1.NSIDE; nside2 = mask2.NSIDE

        # Get good (non-empty, seen) HEAPixels in masks
        good_map_hpix1 = mask1.all_nside_hpix[mask1.seen]
        good_map_hpix2 = mask2.all_nside_hpix[mask2.seen]

        ra_icrs = coords.icrs.ra.deg
        dec_icrs = coords.icrs.dec.deg

        if (mask1.coordframe == 'galactic') | (mask2.coordframe == 'galactic'):
            # Gotta transform...
            galactic_coords = coords.galactic
            ra_galactic = galactic_coords.l.deg
            dec_galactic = galactic_coords.b.deg

        if mask1.coordframe == 'galactic':
            ipix1 = hp.ang2pix(nside1, ra_galactic, dec_galactic,
                        lonlat=True, nest=False)
        else:
            ipix1 = hp.ang2pix(nside1, ra_icrs, dec_icrs,
                        lonlat=True, nest=False)

        if mask2.coordframe == 'galactic':
            ipix2 = hp.ang2pix(nside2, ra_galactic, dec_galactic,
                        lonlat=True, nest=False)
        else:
            ipix2 = hp.ang2pix(nside2, ra_icrs, dec_icrs,
                        lonlat=True, nest=False)

        # Identify whether or not a coordinate lies within the respective
        # masks' seen HEALPixels. seen1/2 are boolean arrays.
        seen1 = np.in1d(ipix1, good_map_hpix1)
        seen2 = np.in1d(ipix2, good_map_hpix2)

        '''
        This returns a boolean array with True if a coordinate is seen in both
        masks, False if its seen in one or neither. NB: gal_ind[seen1 == seen2]
        can also work, but it is risky (what if galaxy is False in both?)
        '''
        overlap = (seen1 == True) & (seen2 == True)

        # Create object index array and return good indices
        gal_ind = np.arange(len(coords))
        good_gals = gal_ind[overlap]

        # Return indices and coodinates of galaxies in both masks
        return good_gals
