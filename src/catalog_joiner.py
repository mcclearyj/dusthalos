import os
from astropy.table import Table, join
import numpy as np
import pdb
from astropy.io import fits
from datetime import datetime

from . import utils, cat_utils

class CatalogJoiner:
    '''
    Helper class to handle the matching of two Catalog objects using either a
    keyword or WCS coordinates. Anticipated usage for CatalogJoiner is being
    instantiated within an instance of Catalog, but CatalogJoiner can also be
    used as a free-standing object.

    Inputs
        config: configuration file specifying paths to each catalog file
        cat1: Catalog instance holding the first catalog
        cat2: Catalog instance holding the second catalog
    '''

    def __init__(self, config, cat1=None, cat2=None):
        self.config = config
        self.cat1 = cat1
        self.cat2 = cat2

    def _cat_loader(self):
        '''
        Placeholder for a method that loads in cat1 if it is a string,
        else does nothing.
        '''

        pass

    def _cat_checker(self):
        '''
        Helper method to convert cat to a Class if it is a dict, then
        checks required arguments using cat_utils.cat_checker() method
        '''

        # A dict is fine too
        try:
            for checkcat in [self.cat1, self.cat2]:
                if (type(checkcat) == dict):
                    checkcat = utils.AttrDict(cat)

                cat_utils.cat_checker(self.cat1)
                cat_utils.cat_checker(self.cat2)

        except AttributeError as ae:
            print('Supplied "cat" is missing required attribute: ', ae)


    def _join_cats_kw(self, overwrite=False):
        '''
        This is a little special-purpose because most catalogs will have photo-z
        and magnitudes in one catalog, but maybe useful anyway.

        Key can be a list!
        '''

        match_kw = self.config['match_kw']
        c1 = self.cat1; c2 = self.cat2
        # Define an output catalog name based on nicknames

        # Join catalogs
        joined_cat = join(c1.data, c2.data,
                        join_type='inner',keys=match_kw,
                        table_names=[c1.tabname, c2.tabname])

        Nmin = min([len(c1.data), len(c2.data)])
        Njoined = len(joined_cat)

        if (Njoined != Nmin):
            print('\n Warning: some galaxies lost during join of photoz ' +
                             'and photometry catalogs;\n' +
                             f'\n\tlen({c1.tabname}) = {len(c1.data)}' +
                             f'\n\tlen({c2.tabname}) = {len(c2.data)}' +
                             f'\n\tlen(joined) = {Njoined}' +
                             '\n'
                             )

        print(f' {Njoined} {c1.tabname} objects joined to {c2.tabname} objects')

        return joined_cat

    def _join_cats_coord(self, cat_dict, overwrite=False):
        '''
        Join catalogs by match in celestial coordinates using
        util.catalog_matcher function.

        Right now, RA and Dec keywords need to be passed explicitly, but could
        read from a config file at some point.

        Note: cat1=None may seem odd, but since we are running this within the
              Catalog() class, can use self.catalog as a default.
        '''

        # Get RA/Dec unit and keys for catalogs 1 and 2;
        # also helps with error checking
        try:
            ra_tag1 = self.cat1.config['ra_tag']
            dec_tag1 = self.cat1.config['dec_tag']
            dec_tag2 = self.config['match']['dec_tag']
            ra_tag2 = self.config['match']['ra_tag']

        except KeyError as ke:
            print('Missing key in config file: ', ke)

        # These are line-by-line matched and can be stacked.
        jcat1, jcat2 = utils.match_coords(self.cat1, self.cat2,
                                ratag1=ra_tag1, dectag1=dec_tag1,
                                ratag2=ra_tag2, dectag2=dec_tag2, radius=0.5)

        # Stack the catalogs
        joined_cat = hstack(jcat1, jcat2, \
                            table_names=[self.cat1.tabname, self.cat2.tabname])

        return joined_cat

    def match_catalogs(self, overwrite=False):
        '''
        Test out whether catalogs are being joined by keyword or RA/Dec,
        call appropriate method for joining
        '''

        # Just for legibility
        config = self.config

        self._cat_checker()

        if (self.config['match_type'] == 'kw'):
            jc = self._join_cats_kw()
            join_type = f'keyword = {self.config["match_kw"]}'

        else:
            jc = self._join_cats_coord()
            join_type = 'coordinates'

        # Add some metadata
        jc.meta['comments'].append(f'{self.cat1.tabname} joined with ' + \
                                    f'{self.cat2.tabname} on' +
                                    f'{datetime.now():%D %H:%M:%S}')
        jc.meta['comments'].append(f'join_type: {join_type}')


        return jc
