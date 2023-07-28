import sys, os
from astropy.table import Table, vstack, hstack
import healpy as hp
import numpy as np
import pdb
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord


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

    def __init__(self, maskfile=None, partial=False, coordframe=None):
        # Do some sanity checking
        if maskfile is None:
            raise ValueError('Missing required parameter "maskfile"')
        if coordframe is None:
            print('Assuming mask is based on equatorial coordinates ("icrs")')
            coordframe = 'icrs'

        # Define attributes
        self.maskfile = maskfile # relative or absolute path
        self.partial = partial
        self.NSIDE = -1
        self.mask = None
        self.mask_header = {}
        self.seen = []
        self.coordframe = coordframe

        # Load mask
        self._load_mask()

    def _load_mask(self):
        '''
        Load up the mask file with a way to handle partial masks
        '''

        try:
            mask, h = hp.read_map(self.maskfile, h=True, partial=self.partial)
        except FileNotFoundError as fnf:
            print(fnf)
        except ValueError as e:
            print(e)
            print("hpMask is probably partial, setting partial=True")
            self.partial = True
            mask, h = hp.read_map(self.maskfile, h=True, partial=self.partial)

        # Set class attributes
        self.mask = mask
        self.mask_header = dict(h)
        self.NSIDE = self.mask_header['NSIDE']
        self.seen = (mask > 0) & (mask != hp.UNSEEN)

        return

    @staticmethod
    def get_overlapping_masks(coords, mask1, mask2):
        '''
        Return coordinates and indices of objects that lie in the overlap area
        of two HEALPix masks. Coords should be an instance of SkyCoords, and
        mask1 and mask2 should be instances of HpMask.
        '''

        # Do some sanity checking
        if type(coords) is not SkyCoord:
            raise TypeError('Supplied "coords" is not an instance of ' + \
                            'astropy.coordinates.SkyCoord')
        if (type(mask1) is not Mask):
            raise TypeError('Supplied "mask1" is not an instance of hpMask')
        if (type(mask2) is not Mask):
            raise TypeError('Supplied "mask2" is not an instance of hpMask')

        nside1 = mask1.NSIDE; nside2 = mask2.NSIDE

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

        # Initialize lists that will hold coordinates and indices of galaxies
        selected_coords = []; good_gals = []

        # Loop
        for j, coord, pix1, pix2 in zip(range(len(coords), coords, ipix1, ipix2):
            if (mask1.mask[pix1]>0) and (mask2.mask[pix2]>0):
                good_gals.append(i); selected_coords.append(coord)

        return good_gals, SkyCoord(selected_coords)
