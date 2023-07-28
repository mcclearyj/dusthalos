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

    def __init__(self, catname, datadir='./', mask=None, config=None):
        # Define attributes
        self.catname = catname
        self.datadir = datadir
        self.mask = mask # should be instance of hpMask
        self.config = config # We may want to define RA tags, DEC tags at some point
        self.catalog = None
        self.joined = None

        self._load_catalog()
        return

    def _load_catalog():
        try:
            catalog = Table.read(os.path.join(self.datadir, self.catname))
        except FileNotFoundError as fnf:
            print(fnf)

        self.catalog = catalog
        return

    @staticmethod
    def join_cats(self, cat1, name1, cat2, name2, key='coadd_object_id', outname='joined.fits', overwrite=False):
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

class FgRandom(Catalog):
    pass
