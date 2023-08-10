import os
from astropy.table import Table, join
import numpy as np
import pdb
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord

import .utils


class Catalog:
    '''
    Hold all catalog information, including joining two different ones and
    identifying overlapping regions based on survey mask.
    '''

    def __init__(self, catfile, path='./', tabname='tab', mask=None, config=None):
        # Define attributes
        self.catfile = catfile
        self.mask = mask # should be instance of hpMask
        self.tabname = tabname #
        self.config = config # We may want to define RA tags, DEC tags at some point
        self.catalog = None
        self.joined_cat = None # In case one is joining an extra catalog

        try:
            catalog = Table.read(os.path.join(self.datadir, self.catfile))
            self.catalog = catalog
        except FileNotFoundError as fnf:
            print(fnf)


    def _cat_checker(self, cat1, cat2, name1, name2):
        '''
        Check that catalogs cat1 and cat2 are either instances of Table or can
        be read in from file. This is an error-handling step useful in joining
        catalogs.
        '''

        # Get catalog 1
        if (cat1 == None):
            cat1 = self.catalog
        elif (type(cat1) == str):
            try:
                cat1 = Table.read(cat1)
            except FileNotFoundError as fnf:
                print(fnf)
        elif (type(cat1) == Table):
            # Don't worry about it
            pass
        else:
            print('cat1 not in a recognized format: str, Catalog, Table')

        # Get catalog2
        if (cat2 == None):
            raise Exception('Parameter "cat2" is None; cannot proceed')
        elif(type(cat2) == Catalog):
            cat2 = cat2.catalog
        elif (type(cat2) == str):
            try:
                cat2 = Table.read(cat2)
            except FileNotFoundError as fnf:
                print(fnf)
        elif (type(cat2) == Table):
            # Don't worry about it
            pass
        else:
            print('cat2 not in a recognized format: str, Catalog, Table')

        if (name1 == None) and (self.tabname == None):
            name1 = 'tab1'
        else:
            name1 = self.tabname

        if (name2 == None):
            name2 = 'tab2'

        return cat1, cat2, name1, name2

    def join_cats_kw(self, cat1=None, name1=None, cat2=None, name2=None,
                        key='coadd_object_id', outname='joined.fits',
                        overwrite=False):
        '''
        This is a little special-purpose because most catalogs will have photo-z
        and magnitudes in one catalog, but maybe useful anyway.

        Key can be a list!
        '''

        # Read in catalogs as needed
        cat1, cat2, name1, name2 = self._cat_checker(cat1, name1, cat2, name2)

        # Join catalogs
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
        joined_cat.write(outname, format='fits', overwrite=overwrite)

        return joined_cat

    def join_cats_coord(self, cat1=None, name1=None, cat2=None, name2=None,
                            ratag1='ra', dectag1='dec', ratag2='ra',
                            dectag2='dec', outname='joined.fits',
                            overwrite=False):
        '''
        Join catalogs by match in celestial coordinates using
        util.catalog_matcher function.

        Right now, RA and Dec keywords need to be passed explicitly, but could be
        read from a config file at some point.

        Note: cat1=None and cat1=None may seem odd, but since we are running this
        within Catalog, it is not always necessary to define additional catalog.
        '''

        # Read in catalogs as needed
        cat1, cat2, name1, name2 = self._cat_checker(cat1, name1, cat2, name2)

        # These are line-by-line matched and can be stacked.
        jcat1, jcat2 = utils.match_coords(cat1, cat2, ratag1=ratag1, dectag1=dectag1,
                                            ratag2=ratag2, dectag2=dectag2,
                                            radius=0.5)
        # Stack the catalogs
        joined_cat = hstack(jcat1, jcat2, table_names=[name1, name2])

        # Save to file, why not
        joined_cat.write(outname, format='fits', overwrite=overwrite)

        return joined_cat


class FgRandom:
    '''
    Subtype of Catalog that just random galaxy catalogs
    '''
    def __init__(self, mask1, config=None):
        pass
