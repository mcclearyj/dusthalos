import numpy as np
import pdb
from src.extinction_model import ExtinctionModel

class ReddeningCalculator(ExtinctionModel):

    def __init__(self, data, redcalc_config=None):
        """
        Do reddening calculation. Though it can be run independently,
        this class is intended to be called by instance of Correlator() class.

        Inherits ExtinctionModel, which holds central wavelengths and instances
        of dust models from the Python 'extinction' library.

        Attributes
            redcalc_config: dust param config; generally suppled as a
                            subset of the top-level run_config in Correlator()
            data: Array-like object with magnitudes, redshifts
            mle: Maximum likelihood estimator for excess reddening
            mle_var: Variance of said MLE
            nbins: Number of redshift bins for Av calculation
            z_key: Column name for self.data redshift info
            use_bin_numbers: Use DES redMaGiC redshift bin numbers? Prob. not.
        """

        self.data = data #
        self.mle = [] #
        self.mle_var = [] # Variance of said MLE
        self.nbins = None # Number of redshift bins for Av calculation
        self.z_key = None # Column name for self.data redshift info

        # dust config should be in here now too
        self.load_config(redcalc_config)

        super().__init__(dmconfig = redcalc_config['dust_params'])

    def load_config(self, redcalc_config):
        """
        redshift config should have all the parameters keywords
        """
        # Resonable redshift defaults, not sure about dust params
        #default_dust_pars = {'model': 'calzetti00', 'R': 3.1, 'A_V': 1,
        #                    'wavelengths': [4808.49, 6417.65, 7814.58, 9168.85]}
        # No, it should fail loudly.

        # Set sensible defaults
        base_config = {
            'nbins': 7, 'z_key': 'z', 'use_bin_numbers': False,
            'dust_params': None
        }

        # Overwrite defaults, append non-standard keys b/c who cares.
        if redcalc_config != None:
            config_keys = redcalc_config.keys()
            for key in config_keys:
                if key not in base_config.keys():
                    print(f'Warning: "{key}" is not a standard ',
                            'ReddeningCalculator config key:')
                    print(f'\t {base_config.keys()}')
                base_config[key] = redcalc_config[key]

        self.redcalc_config = base_config
        self.dust_params = base_config['dust_params']

        print(f'\n ReddeningCalculator() configuration values:')
        print(f'\t {base_config}')

    def preprocess_catalog(self):
        """
        If there are NaNs or masked values or whatever, clean them out of the
        catalog and be noisy about it.
        """
        # Placeholder bandholder
        bands = self.dust_params['band_names']

        # Start by assuming all gals in bands are good <3
        catlen = len(self.data)
        wg = np.full(catlen, True)

        # Pick out only the good galaxies!
        for band in bands:
            band_mag = np.ma.getdata(self.data[band])
            band_bool = (band_mag > -9999) & \
                        (band_mag != np.nan) & \
                        (band_mag < 28)
            wg *= band_bool

        # percent of galaxies that failed to pass selections
        pfail = 100-(np.count_nonzero(wg) / catlen * 100)

        # How many failed?
        print(f'ReddeningCalc: {np.count_nonzero(wg)}/{catlen} galaxies',
              f'({100-pfail:.1f}%) have good photometry')
        print(f'ReddeningCalc: removing {pfail:.1f}% of galaxies from data')
        print('')

        self.data = self.data[wg]
        self.good_indices = wg

    def optimal_estimator(self, catalog):
        """
        Compute the optimal estimator for excess reddening based on the
        input galaxy catalog.

        TO DO: Move error-checking to earlier in the code?
        """

        # Extract galaxy magnitudes from the catalog
        band_names = self.dust_params['band_names']
        data_list = []
        for band in band_names:
            try:
                # Append each magnitude column to the data object
                data_list.append(np.ma.getdata(catalog[band]))
            except KeyError as ke:
                # Oh no! Exit, as there are probably other problems
                print(f"ReddeningCalcuator: no bandpass named '{band}' found!")
                raise ke

        # Stack the extracted magnitudes into a 2D array
        data = np.vstack(data_list)

        # Compute the covariance matrix of the galaxy magnitudes
        covar = np.cov(data)

        # Initialize delta based on the dust extinction model (self.dmdp)
        delta = (np.zeros_like(data).T + self.dmdp)

        # Compute the inverse covariance matrix
        Cinv = np.linalg.inv(covar)

        # De-mean the galaxy magnitudes in each redshift bin
        colors = (data.T - np.median(data, axis=1)).T

        # Compute the optimal estimator for excess reddening of galaxies
        # in this redshift bin using maximum likelihood
        est = np.sum(delta.T * np.dot(Cinv, colors), axis=0) / \
                        np.sqrt(np.dot(self.dmdp, np.dot(Cinv, self.dmdp)))

        # Compute the Cramer-Rao bound for the optimal estimator
        wt = 1./np.sqrt(np.dot(self.dmdp, np.dot(Cinv, self.dmdp)))

        return est, wt

    def run(self):
        """
        zmin: minimum redshift of galaxies to consider (default None)
        nbins: number of bins for color estimation
        """

        print(f'Removing catalog entries with NaN & sentinel values')
        self.preprocess_catalog()

        print(f'Make dust extinction model')
        self.get_dust_model()

        print(f'Beginning background catalog reddening calculation...')

        mle = np.zeros(len(self.data))
        mle_var = np.zeros(len(self.data))

        if self.redcalc_config['use_bin_numbers'] == True:
            print("using bin numbers for redshifts")
            zbin_col = self.data['bin_number']
        else:
            # Make a "bin number" on the fly!
            print("Making redshift bins")
            z_hist = np.histogram(self.data[self.redcalc_config['z_key']],
                                    bins=self.redcalc_config['nbins'])
            zbin_col = np.digitize(self.data[self.redcalc_config['z_key']],
                                    bins=z_hist[1], right=False)

        bin_numbers = np.unique(zbin_col)

        try:
            # Loop over all redshift bins specified in bin_numbers
            for zb in bin_numbers:
                # Create a boolean mask to select only galaxies in
                # the current redshift bin
                slice = (zbin_col == zb)

                # Compute the optimal estimator and Cramer-Rao bound for
                # galaxies in the current redshift bin
                this_est, this_wt = self.optimal_estimator(self.data[slice])

                # Assign the computed optimal estimator and Cramer-Rao
                # bound to the output arrays
                mle[slice] = this_est
                mle_var[slice] = this_wt

        except np.linalg.LinAlgError:
            # If a Linear Algebra Error occurs, it's likely due to having too
            # few galaxies in the bin
            error = f'Too few galaxies in bin {zb}: {np.count_nonzero(slice)}'
            mle[slice] = np.nan
            mle_var[slice] = np.nan

        self.mle = mle
        self.mle_var = mle_var

    def make_treecorr_obj(self):
        """
        Placeholder for function that can return a treecorr object with k, w if desired!
        """
        self.treecorr_cat = treecorr.Catalog(ra=self.coords.ra.deg,
                            dec=self.coords.dec.deg, ra_units='deg',
                            dec_units='deg', k=redcalc.mle,
                            w=redcalc.mle_var)
