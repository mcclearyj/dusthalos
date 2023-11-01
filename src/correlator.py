import numpy as np
import os
import pdb
import time
import treecorr
#matplotlib.use('Agg')
from astropy.io import fits
from astropy.table import Table

# Local imports
#from src.dust_calculator import DustCalculator
import src.utils as utils
from src.catalog import Catalog
from src.cat_utils import cat_config_checker
from src.reddening_calc import ReddeningCalculator

class Correlator:
    '''
    Load & store all the Treecorr.Catalogs to be used. Also run the
    reddening calculations for the background objects. Actually, wait, the
    Catalog object have everything I need, including utilities for grabbing &
    transforming coordinates. Just use those!

    correl config should hold... what, locations of cat configs?
    Certainly where to store them. First test should just be can we load the
    configuration files in the first place
    '''

    def __init__(self, correl_config=None, ctype=None):
        self.correl_config = correl_config
        self.Catalog = None
        self.treecorrCatalog = None
        self.ctype = ctype

    def _read_cat_file(self, path=None):
        '''
        Read in FITS catalog as a Table() object
        '''
        if path is None:
            # load catalog
            pass
        else:
            pass

    def _load(self):
        '''
        Parameter
            cat_config: a Catalog-format config file and create Catalog object
        Returns:
            treecorr_cat: treecorr.Catalog object
        '''
        this_cat = Catalog(config=self.cat_config)
        coords = this_cat.grab_coords().icrs # Returns SkyCoord instance
        self.coords = coords
        tc_cat = treecorr.Catalog(ra=coords.ra.deg, dec=coords.dec.deg,
                                ra_units='deg',dec_units='deg')
        return this_cat, tc_cat

    def load(self):
        '''
        Parameter
            cat_config: a Catalog-format config file and create Catalog object
        Returns:
            treecorr_cat: standard treecorr.Catalog object, no extra frink
        '''
        # Can't run in the catalog type or whatever isn't in correl_config
        ctype = self.ctype
        assert self.ctype in self.correl_config.keys()

        # Error checking
        cat_config_path = self.correl_config[ctype]['cat_config']
        self.cat_config = utils.read_yaml(cat_config_path)
        cat_config_checker(self.cat_config)

        # Set Catalog and treecorrCatalog objects
        self.Catalog, self.treecorrCatalog = self._load()

    def do_reddening(self):
        '''
        Call ReddeningCalculator, run it
        '''
        ctype = self.ctype

        # Additional checker: make sure we have redshift
        if 'redshifts' not in self.correl_config[ctype].keys():
            raise KeyError('correl_config missing parameter group "redshifts"')
        else:
             rc_config = self.correl_config[ctype]['redshifts']

        if 'z_tag' not in self.cat_config.keys():
            raise KeyError('catalog_config missing parameter "z_tag" ')
        else:
            rc_config['z_tag'] = self.cat_config['z_tag']

        # Grab dust parameters
        dust_model_config = self.correl_config['dust_params']
        rc_config['dust_params'] = dust_model_config

        # Instantiate ReddeningCalculator
        reddening_calc = ReddeningCalculator(self.Catalog.data,
                                        redcalc_config=rc_config)
        reddening_calc.run()

        # Grab indices with clean photometry
        wg = reddening_calc.good_indices

        # I can't figure out a good way to organize treecorr.Catalog creation
        updated_treecorr_catalog = treecorr.Catalog(ra=self.coords[wg].ra.deg,
                                    dec=self.coords[wg].dec.deg, ra_units='deg',
                                    dec_units='deg', k=reddening_calc.mle,
                                    w=reddening_calc.mle_var)

        updated_treecorr_catalog.redshift = \
                    self.Catalog.data[wg][self.cat_config['z_tag']]

        self.treecorrCatalog = updated_treecorr_catalog

    def write_to_file(self, outname=None):
        '''
        Save Treecorr catalog to file
        '''
        if outname == None:
            outname = os.path.join(self.correl_config['output_path'],
                    self.correl_config['output_basename']+'_treecorrcat.fits')
        self.treecorrCatalog.write(outname)


    def load_all_cats(self):
        '''
        Something of a utility function to step through and populate all the
        treecorr catalogs with something
        '''

        pass
