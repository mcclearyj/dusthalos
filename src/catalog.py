import os
from astropy.table import Table, join
import numpy as np
import pdb
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from datetime import datetime

# Local imports
from . import utils, cat_utils
from .hpmask import HpMask
from .catalog_joiner import CatalogJoiner

class Catalog():
    '''
    Hold all catalog information, including joining two different ones and
    identifying overlapping regions based on survey mask. Inherits CatalogJoiner
    and CatalogMasker classes.
    -----------------
    Inputs/Attributes
            filename: (base) name of input catalog

            path: path to input catalog [default='./']

            tabname: nickname for table if matching to another catalog and
                     the two tables have identical column names [default='tab']

            mask: mask file (path name or HpMask instance) for this catalog

            config: sets attributes of catalog class

            memmap: "secret" argument that enables memmap for large fits
                    files -- should really be made a configurable parameter.
    '''

    def __init__(self, config=None, filepath=None, data=None,
                 tabname='tab', mask=None, vb=True, memmap=False):
        self.config = config # configuration file specifying catalog parameters
        self.filepath = filepath # absolute or relative path to catalog file
        self.mask = mask # should be instance of hpMask
        self.tabname = tabname # can be convenient if matching
        self.data = data # Table data for catalog
        self.vb = vb # Verbosity

        self._load_config()

        # If filepath to catalog has been set, load it
        if self.filepath is not None:
            self._load_cat(self.filepath, memmap=memmap)

    def _load_config(self):
        '''
        If it exists, load the config
        '''
        if self.config is not None:
            self.filepath = os.path.join(self.config['path'],
                                self.config['filename'])
            if 'tabname' in self.config.keys():
                self.tabname = self.config['tabname']
        else:
            pass

    def _load_cat(self, input, memmap):
        '''
        Utility to load in a catalog and set attributes of Catalog() object
        then populates self.data and self.tabname
        Input
            input:  one of filepath or astropy Table
            memmap: whether or not to use memmap (for large tables)
        '''

        # Sanity check
        if type(input) not in [str, Table]:
            raise TypeError(f'input of type {type(input)} not in [str, Table]')

        # Load catalog data
        if type(input) == str:
            try:
                if memmap == True:
                    catalog_data = Table.read(input, format='fits', memmap=True)
                else:
                    catalog_data = Table.read(input)
                print(f'Loaded catalog data from {input}')

            except FileNotFoundError as fnf:
                print('Invalid catalog path supplied: ', fnf)
        else:
            catalog_data = input
            print(f'Setting catalog data to input Table')

        # Set attributes
        self.data = catalog_data

        # Also intialize a meta file
        self.data.meta = \
            {'comments': [f'Original catalog length = {len(self.data)}']}

    def _write_outfile(self, prefix, overwrite=True):
        # Defile output file name & save
        outcat_name = f'{prefix}{os.path.basename(self.filepath)}'
        outcat_path = os.path.join(self.config['output_path'], outcat_name)
        self.data.write(outcat_path, format='fits', overwrite=overwrite)
        if self.vb == True:
            print(f'\tSaved catalog to {outcat_path}')

    def grab_coords(self):
        '''
        Utility function to return a SkyCoord instance fromt the coordinates
        of galaxies in self.data
        '''
        # What format are the catalog's coordinates in?
        if self.config['coord_units'] == 'deg':
            unit = u.deg
        else:
            unit = u.rad

        # Create and return SkyCoord instance
        coords = SkyCoord(self.data[self.config['ra_key']].data*unit,
                            self.data[self.config['dec_key']].data*unit,
                            frame=self.config['coordframe'])

        return coords

    def create_mask_from_config(self, mask_config=None):
        '''
        Method to create an HpMask object from the Catalog configuration file

        Input
            mask_config: specifies mask basename, absolute/relative path,
                         and coordframe (icrs or galactic)
        Returns
            self.mask: HpMask object for this catalog
        '''

        mask_filepath = os.path.join(mask_config['path'],
                                     mask_config['filename'])

        kw_args = {'coordframe': mask_config['coordframe']}

        if 'partial' in mask_config.keys():
            kw_args['partial'] = mask_config['partial']

        self.mask = HpMask(mask_filepath, **kw_args)


    def apply_mask(self, overwrite=True):
        '''
        Method to apply a HEAPix mask using the HpMask.apply_mask method!
        Trims self.data so that only objects whose coordinates fall inside
        one of the allowed HEALPixels remain.
        '''

        # Get current catalog length
        curr_len = len(self.data)

        # Make sure mask exists first!
        if (type(self.mask) is not HpMask):
            raise TypeError("Can't apply a mask of type: {type(self.mask)}")

        # Get catalog's coordinates
        coords = self.grab_coords()

        # Obtain indices of seen galaxies using HpMask.apply_mask()
        good_gals = self.mask.apply_mask(coords)

        # Set data to masked data
        self.data = self.data[good_gals]

        # update meta information
        self.data.meta['comments'].append(
                    f'Applied HEALPix mask {self.mask.filepath} on ' + \
                    f' date {datetime.now():%D %H:%M:%S}')
        self.data.meta['comments'].append(
                    f'masked length: {len(self.data)}')

        print(' Catalog.apply_mask(): ')
        print('\tSet self.data to healpix-masked catalog')

        # Write outfile
        self._write_outfile(prefix='Masked_')

    def apply_overlapping_masks(self, mask1=None, mask2=None):
        '''
        Runs HpMask.find_overlapping_masks() with mask1 = self.mask
        and mask2 passed as an argument.
        '''

        # Hold original, unmasked catalog data
        curr_len = len(self.data)

        # Sensible default
        if mask1 == None:
            mask1 = self.mask

        # Grab catalog coords
        coords = self.grab_coords()

        # Run HpMask.find_overlapping_masks()
        good_gals = HpMask.apply_overlapping_masks(coords, mask1, mask2)

        # Update self.data to reflect masking and update meta information
        self.data = self.data[good_gals]

        self.data.meta['comments'].append(\
                    f'Applied HEALPix masks {mask1.filepath} and ' + \
                    f'{mask2.filepath} on date {datetime.now():%D %H:%M:%S}')
        self.data.meta['comments'].append(\
                    f'masked length: {len(self.data)}')

        # Print a reassuring message
        print(' Catalog.apply_overlapping_masks(): ')
        print('\tSet self.data to objects at intersection of mask1, mask2')

        # Write outfile
        self._write_outfile(prefix='DoubleMasked_')

    def match_to_catalog(self, catalog2, overwrite=False):
        '''
        Calls CatalogJoiner to joins "catalog1" (self) to the catalog2
        passed as argument. Catalog2 should have a Catalog-type format with data
        and tabname attributes.
        '''

        # CatalogJoiner does all the checking and stuff
        joiner = CatalogJoiner(self.config['match'], cat1=self, cat2=catalog2)

        # Return the matched/joined
        jc = joiner.match_catalogs()

        # Define an output catalog name based on nicknames and write to file
        output_name = \
            f'{self.tabname}_{catalog2.tabname}_JOINED_catalog.fits'
        output_file = os.path.join(self.config['output_path'], output_name)
        jc.write(output_file, format='fits', overwrite=overwrite)

        # Overwrite existing catalog with joined catalog
        self.data = jc

        if self.vb is True:
            print('\n Setting self.data to joined catalog \n')
