import sys, os
from astropy.table import Table, vstack, hstack
import healpy as hp
import numpy as np
import pdb
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt; plt.ion()

import .utils

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


class BkgCat():
    '''
    Placeholder
    '''

    def __init__(self, maskfile=None, photo_catalog_name=None,
                redshift_catalog_name=None):
        # Define attributes
        self.photo_catalog_name = photo_catalog_name
        self.redshift_catalog_name = redshift_catalog_name
        self.output_catalog_name = outcat
        self.photo_catalog = None
        self.redshift_catalog = None
        self.joined_catalogs = None
        self.z_min = None

        return


    def do_mof_selections(gold catalog):
        # gold_catalog should be a table object
        assert(rmz_joined)
        keep = (rmz_joined['mof_flags']==0) & \
        (rmz_joined['mof_cm_flux_g']>0) & (rmz_joined['mof_cm_flux_r']>0) & \
        (rmz_joined['mof_cm_flux_i']>0) & (rmz_joined['mof_cm_flux_z']>0) & \
        (rmz_joined['zredmagic'] > zmin)

        return


class FgRandoms(hpMask):
    '''
    *** DOES NOTHING RIGHT NOW. SEE get_fg_random.py FOR HOW THEY ARE MADE ***

    FgRandoms class represents a set of foreground random points.

    This class inherits from the hpMaskMixin class to provide methods
    for handling masks and regions of interest.

    Parameters:
        masks (list): List of masks to apply to the random points.
                      Default is an empty list.

        outname (str): Name of the output file for the random points.
                       Default is 'foreground_random.fits'.

        overwrite (bool): Flag to indicate whether to overwrite an existing
                          output file. Default is True.

        outdir (str): Output directory for the file. Default is './'.
    '''

    def __init__(self, masks=None,
                outname='foreground_random.fits',
                overwrite=True, outdir='./'):

        self.masks = masks
        self.outname = outname
        self.outdir = outdir
        self.overwrite = overwrite


    def _make_fg_randoms_old(self, maskfile,
                            nrand=1e6, nside=4096,
                            partial=True, nest=False):
        # Make randoms on the sphere. Deprecated.
        maskfile = self.maskfiles
        print('loading hmap, making randoms')

        hmap, hd = hp.read_map(maskfile, partial=partial, nest=nest, h=True)
        # But how many? Try to get approximately nrand, if possible.
        fcover = np.sum(hmap > 0)*1./hmap.size
        ndraw = np.ceil(nrand/fcover*1.2).astype(int)

        ran1, ran2 = np.random.random(2*ndraw).reshape(2, -1)
        ra  = 2. * np.pi * (ran1 - 0.5) * 180./np.pi
        dec = np.arcsin(2. * (ran2 - 0.5)) * 180./np.pi

        hpInd = hpRaDecToHEALPixel(ra,dec,nside=nside,nest=nest)
        keep  = hmap != hp.UNSEEN

        use = np.random.rand(ra.size) < hmap[hpInd]

        ra  = ra[use]
        dec = dec[use]
        ra = ra[:int(nrand)]
        dec = dec[:int(nrand)]
        # covert ra from [-180,180 )  to [0,360)
        ra = (ra + 360) % 360

        fg_out = Table([ra, dec], names=['ra', 'dec'])
        fg_out.write('wiseSCOS_random_catalog.fits', overwrite=True)
        # Need to double-check this RA/Dec conversion

        print('making fg random treecorr catalog\n')

        fgRan = treecorr.Catalog(ra=ra, dec=dec,
                ra_units='deg',dec_units='deg')

        return fgRan

class Catalog:
    '''
    Hold all catalog information, including joining two different ones and
    identifying overlapping regions based on survey mask.
    '''

    def __init__(cat_name=None):
        # Define attributes
        self.cat_name = None
        self.joined = None
        self.thing = None



    def join_cats(cat1, name1, cat2, name2, key='coadd_object_id', outname='joined.fits', overwrite=False):
        '''
        This is a little special-purpose because I don't anticipate needing to join
        photo-z and magnitude catalogs a lot, but maybe useful anyway.

        Key can be a list!

        # rename a few cols for joining purposes
        colmap = {
            'ALPHAWIN_J2000': 'ra',
            'DELTAWIN_J2000': 'dec',
            'NUMBER': 'id'
            }
        for old, new in colmap.items():
            self.det_cat.rename_column(old, new)
        '''

        joined = join(cat1, cat2, join_type='inner',
                    keys=[key], table_names=[name1, name2])

        Nmin = min([len(cat1), len(cat2)])
        Nobjs = len(joined)

        if (Nobjs != Nmin):
            print('There was an error while joining the photoz and ' +
                             'photometry catalogs;' +
                             f'\nlen({name1})={len(cat1)}' +
                             f'\nlen({name2})={len(cat2)}' +
                             f'\nlen(joined)={Nobjs}'
                             )

        print(f'{Nobjs} {name1} objects joined to {name2} objects')
        joined.write(outname, format='fits', overwrite=overwrite)

        return joined

    def get_overlapping_masks(coords, mask1, mask2):
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

        selected_coords = []; good_gals = []
        for j, coord, pix1, pix2 in zip(range(len(coords), coords, ipix1, ipix2):
            if (mask1.mask[pix1]>0) and (mask2.mask[pix2]>0):
                good_gals.append(i); selected_coords.append(coord)

        return good_gals, SkyCoord(selected_coords)

def match_DES_catalogs(datapath, phot_file, rmz_file, zmin=0.15):
    '''
    Utility function for creating a DES background catalog that also includes
    MOF photometry. Joining is based on COADD ids;
    ALSO DON'T FORGET TO MASK BOTH!
    '''
    # Doing a join would be more solid
    bg_phot = Table.read(phot_file)
    bg_redshift = Table.read(rmz_file)
    outname = 'joined.fits'
    print 'phot and z files read in'

    # Match these
    joined = join_cats(bg_phot, 'Y3gold', bg_redshift, 'redmagic',
                        outname=outname, overwrite=True)

    print "cats matched"

    # Now for the fun part... mask them


    print "background redMaGiC catalog acquired"

    return redMaGiC
