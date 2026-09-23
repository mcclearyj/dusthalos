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

    Attributes
        correl_config: high-level run configuration file
                       see 'configs/dust_calc_config.yaml' for an example
        cat_config: config file with locations of data catalog, colnames, etc.
        ctype: type of catalog for calculation; should be one of
               'background_catalog', 'background_randoms',
               'foreground_catalog', 'foreground_randoms'.
        Catalog: instance of dust_halos Catalog class that holds mask,
                 catalog data, etc. Left as None for catalogs loaded straight
                 from a saved treecorr catalog.
        treecorrCatalog: instance of treecorr Catalog class, which is returned
                         by reddening_calc(). Contains Av mle, wt_mle
        from_file: True if treecorrCatalog was read from a saved treecorr
                   catalog rather than built from a science catalog. In that
                   case the k values are already in the file, so do_reddening()
                   and write_to_file() are no-ops.
    '''

    # Defaults match the column names written by treecorr.Catalog.write()
    treecorr_cat_defaults = {
        'ra_col': 'ra',
        'dec_col': 'dec',
        'w_col': 'w',
        'k_col': 'k',
        'patch_col': 'patch',
        'ra_units': 'deg',
        'dec_units': 'deg'
    }

    def __init__(self, correl_config=None, ctype=None):
        self.correl_config = correl_config
        self.cat_config = None
        self.Catalog = None
        self.treecorrCatalog = None
        self.ctype = ctype
        self.coords = None
        self.from_file = False

    def _read_cat_file(self, path=None):
        '''
        Read in FITS catalog as a Table() object
        '''
        if path is None:
            # load catalog
            pass
        else:
            pass

    def load(self, treecorr_npatch=None, treecorr_patch_centers=None):
        """
        Populate dusthalos-type Catalog attribute (self.catalog) based on
        config values after checking it is one of background, background_random,
        foreground, foreground_random. Also populates the treecorr-type Catalog
        attribute (self.treecorrCatalog).

        If the run config supplies a 'treecorr_cat' block for this ctype, the
        treecorr catalog is read directly from that file instead (see
        _load_from_file), and no dusthalos Catalog is created. Catalog types
        with and without a 'treecorr_cat' block can be mixed freely in one run.
        """
        # Can't load a dusthalo-catalog type not specified in correl_config
        ctype = self.ctype
        assert self.ctype in self.correl_config.keys()

        # Saved treecorr catalog? Then there's no science catalog to read
        file_config = self.correl_config[ctype].get('treecorr_cat')

        if file_config is not None:
            self.from_file = True
            self.treecorrCatalog = self._load_from_file(
                file_config, treecorr_npatch, treecorr_patch_centers
            )
            return

        # Check configuration values for errors
        cat_config_path = self.correl_config[ctype]['cat_config']
        self.cat_config = utils.read_yaml(cat_config_path)
        cat_config_checker(self.cat_config)

        # If a treecorr_npatch argument isn't passed, look for one in config
        # If there isn't one, returns "None"
        if (treecorr_npatch == None):
            treecorr_npatch = self.correl_config[ctype].get('treecorr_npatch')

        # Set Catalog and treecorrCatalog attributes
        self.Catalog, self.treecorrCatalog = \
            self._load(treecorr_npatch, treecorr_patch_centers)

    def _load(self, treecorr_npatch, treecorr_patch_centers):
        '''
        Parameter
            cat_config: a Catalog-format config file and create Catalog object
        Returns:
            treecorr_cat: treecorr.Catalog object
        '''
        this_cat = Catalog(config=self.cat_config)
        coords = this_cat.grab_coords().icrs # Returns SkyCoord instance
        self.coords = coords

        kwargs = {
            'patch_centers': treecorr_patch_centers,
            'ra_units': 'deg',
            'dec_units': 'deg'
        }

        if treecorr_npatch != None:
            kwargs['npatch'] = int(treecorr_npatch)

        tc_cat = treecorr.Catalog(
            ra=coords.ra.deg, dec=coords.dec.deg, **kwargs
        )

        return this_cat, tc_cat

    def _load_from_file(self, file_config, treecorr_npatch=None,
                            treecorr_patch_centers=None):
        '''
        Build a treecorr.Catalog from an already-saved treecorr catalog, e.g.
        one written by Correlator.write_to_file() and de-meaned by
        DemeanTreecorrCat. The k values are taken from the column named by
        'k_col', which is how one selects e.g. 'demeaned_k' over 'k'.

        Parameters
            file_config: 'treecorr_cat' block of the run config for this ctype.
                         Only 'filename' is required; column names and coord
                         units default to Correlator.treecorr_cat_defaults.
            treecorr_npatch: only used if the file has no patch column
            treecorr_patch_centers: only used if the file has no patch column
        Returns:
            treecorr_cat: treecorr.Catalog object
        '''
        if 'filename' not in file_config.keys():
            raise KeyError(
                f"treecorr_cat config for '{self.ctype}' is missing " + \
                "required key 'filename'"
            )

        fname = file_config['filename']

        if not os.path.isfile(fname):
            raise FileNotFoundError(
                f"treecorr_cat file for '{self.ctype}' not found: {fname}"
            )

        # Fill in any column names & units the config didn't supply
        params = dict(self.treecorr_cat_defaults)
        params.update(
            {k: v for k, v in file_config.items() if k != 'filename'}
        )

        # Which columns does the file actually have? Header only, no data read:
        # these catalogs run to several GB, so let treecorr do the row loading
        with fits.open(fname, memmap=True) as hdul:
            colnames = hdul[1].columns.names

        if params['k_col'] not in colnames:
            raise KeyError(
                f"treecorr_cat file {fname} has no column " + \
                f"'{params['k_col']}'; available columns are {colnames}"
            )

        # Weights and patches are optional; drop them if they aren't in there
        for optional_col in ['w_col', 'patch_col']:
            if params[optional_col] not in colnames:
                print(f"  no '{params[optional_col]}' column in {fname}, " + \
                        f"ignoring {optional_col}")
                params[optional_col] = None

        # A patch column fixes the patch assignment, and treecorr won't accept
        # it alongside npatch/patch_centers. Without one, fall back to the
        # patch arguments so the catalog still shares the run's patch geometry
        if params['patch_col'] is not None:
            if (treecorr_npatch is not None) \
                    or (treecorr_patch_centers is not None):
                print(f"  {self.ctype}: using patch column " + \
                        f"'{params['patch_col']}' from file; ignoring " + \
                        "supplied npatch/patch_centers")
        else:
            params['patch_centers'] = treecorr_patch_centers
            if treecorr_npatch is not None:
                params['npatch'] = int(treecorr_npatch)

        print(f"Loading treecorr catalog for {self.ctype} from {fname}")
        print(f"  using k_col='{params['k_col']}', " + \
                f"patch_col={params['patch_col']}")

        return treecorr.Catalog(file_name=fname, **params)


    def do_reddening(self):
        """ 
        Check dust_params config, call ReddeningCalculator, run it. The various 
        configs do get a little messy, and it might be worth reorganizing them. 
        For now: the 'correl_config' is the overall dust_calc_config, organized
        by 'ctype' (bg random, bg galaxy, fg rand, fg galaxy). 'cat_config' is
        the  catalog_config with all relevant column names and file paths for 
        that particular catalog (redMaGiC, WISE, ...)
        
        """

        ctype = self.ctype

        # A catalog read from file already has its k (and w) values
        if self.from_file == True:
            print(f"{ctype} was loaded from a saved treecorr catalog, " + \
                    "skipping reddening calculation")
            return

        # Additional checker: make sure we have redshift
        if 'redshifts' not in self.correl_config[ctype].keys():
            raise KeyError(f"run config missing parameter set 'redshifts'")
        else:
            rc_config = self.correl_config[ctype]['redshifts']

        # Check redshifts
        if 'z_key' not in self.cat_config.keys():
            raise KeyError("catalog_config missing parameter 'z_key'")
        else:
            rc_config['z_key'] = self.cat_config['z_key']

        # Update bin number parameters if using them in Av calculation
        if ('use_bin_numbers' in rc_config.keys()) == True:
            if 'zbin_key' not in self.cat_config.keys():
                raise KeyError(
                    "catalog_config missing parameter 'zbin_key' required " + \
                    "for using redshift bin numbers in A_V calculation"
                )
            else:
                rc_config['zbin_key'] = self.cat_config['zbin_key']
            
        # Grab dust parameters
        try:
            dust_model_config = self.correl_config['dust_params']
            rc_config['dust_params'] = dust_model_config
        except KeyError as ke:
            print(f"run config missing parameter set 'dust_params'")

        # Instantiate ReddeningCalculator
        reddening_calc = ReddeningCalculator(
            self.Catalog.data, redcalc_config=rc_config
        )
        reddening_calc.run()

        # Grab indices with clean photometry
        wg = reddening_calc.good_indices

        # Grab patch centers from OG catalog
        patch_centers = self.treecorrCatalog.patch_centers

        # I can't figure out a good way to organize treecorr.Catalog creation
        updated_treecorr_catalog = treecorr.Catalog(
            ra=self.coords[wg].ra.deg,
            dec=self.coords[wg].dec.deg, ra_units='deg',
            dec_units='deg', k=reddening_calc.mle,
            w=reddening_calc.mle_var,
            patch_centers = patch_centers
        )

        # Add redshift, too
        updated_treecorr_catalog.redshift = \
                    self.Catalog.data[wg][self.cat_config['z_key']]

        self.treecorrCatalog = updated_treecorr_catalog

    def write_to_file(self, outname=None):
        '''
        Save Treecorr catalog to file
        '''
        # Don't overwrite the file this catalog was just read from
        if self.from_file == True:
            print(f"{self.ctype} was loaded from a saved treecorr catalog, " + \
                    "not writing it back out")
            return

        if outname == None:
            outname = os.path.join(self.correl_config['output_path'],
                    self.correl_config['output_basename']+'_treecorrcat.fits')
        self.treecorrCatalog.write(outname)

