import numpy as np
import pdb

class ReddeningCalculator:

    def __init__(self, data, redshift_config=None):
        '''
        Do reddening calculation
        '''
        self.data = data # Data with magnitudes, redshifts
        self.mle = [] # Maximum likelihood estimator for anomalous reddening
        self.mle_var = [] # Variance of said ML
        self.nbins = None
        self.z_tag = None
        self.use_bin_numbers = False

        self.load_config(redshift_config)

    def load_config(self, redshift_config=None):
        '''
        redshift config should have all the parameters keywords
        '''

        # Resonable defaults
        default_config = {'nbins': 7,
                          'z_tag': 'z',
                          'use_bin_numbers': False,
                          'dust_model': 'ccm89',
                          'bandpasses': None
                          }

        # If you start adding others, they won't do much
        allowed_keys = default_config.keys()

        # Set default
        if redshift_config == None:
            self.config = default_config

        else:
            config_keys = redshift_config.keys()

            for key in config_keys:
                # Overwrite defaults, append non-standard keys b/c who cares.
                if key not in allowed_keys:
                    print(f'"{key}" is not a standard ReddeningCalculator config key:')
                    print(f'\t {allowed_keys}')
                default_config[key] = redshift_config[key]

            self.config = default_config

        print(f'Using ReddeningCalculator configuration values:')
        print(f'\t {default_config}')

    def preprocess_catalog(self):
        '''
        If there are NaNs or masked values or whatever, clean them out of the
        catalog and be noisy about it.
        '''

        # Placeholder bandholder
        bands = ['mof_cm_mag_corrected_g', 'mof_cm_mag_corrected_r',
                'mof_cm_mag_corrected_i', 'mof_cm_mag_corrected_z']

        # All bands are good <3
        catlen = len(self.data)
        wg = np.full(catlen, True)

        # Pick out only the good entries!
        for band in bands:
            band_mag = np.ma.getdata(self.data[band])
            band_bool = (band_mag > -9999) & (band_mag != np.nan)
            wg *= band_bool

        # percent of galaxies that failed
        pfail = 100-(np.count_nonzero(wg) / catlen * 100)

        # How many failed?
        print(f'RedshiftCalc: {np.count_nonzero(wg)}/{catlen} galaxies',
              f'({100-pfail:.1f}%) have good photometry')
        print(f'RedshiftCalc: removing {pfail:.1f}% of galaxies from data')
        print('')

        self.data = self.data[wg]
        self.good_indices = wg

    def optimal_estimator(self, catalog):
        '''
        Make the colors. Catalog... data... too many uses of same thing
        TO DO:
            - make bandpass names configurable
            - be able to select a specific dust model
        '''

        # Fitzpatrick99
        fitz99 = np.array([1.12150099, 0.77164321, 0.57725486, 0.45259124])
        # Calzetti:
        calz = np.array([1.13552323, 0.77956032, 0.55082914, 0.39834168])
        # Cardelli, Clayton & Mathis '89
        ccm89 = np.array([1.12224688, 0.82747095, 0.62680647, 0.47880753])

        dmdp = ccm89

        gmag = np.ma.getdata(catalog['mof_cm_mag_corrected_g'])
        rmag = np.ma.getdata(catalog['mof_cm_mag_corrected_r'])
        imag = np.ma.getdata(catalog['mof_cm_mag_corrected_i'])
        zmag = np.ma.getdata(catalog['mof_cm_mag_corrected_z'])

        data = np.vstack([gmag,rmag,imag,zmag])
        covar = np.cov(data)

        delta = (np.zeros_like(data).T + dmdp)
        Cinv = np.linalg.inv(covar)
        colors = (data.T - np.average(data, axis=1)).T
        est = np.sum(delta.T * np.dot(Cinv,colors),axis=0) / \
                        np.sqrt(np.dot(dmdp, np.dot(Cinv, dmdp)))
        wt = 1./np.sqrt(np.dot(dmdp, np.dot(Cinv, dmdp)))

        return est, wt

    def calc_reddening(self):
        '''
        zmin: minimum redshift of galaxies to consider (default None)
        nbins: number of bins for color estimation
        '''

        print(f'Removing catalog entries with NaN & sentinel values')
        self.preprocess_catalog()

        print(f'Beginning background catalog reddening calculation...')

        nbins = self.config['nbins']
        z_tag = self.config['z_tag']

        mle = np.zeros(len(self.data))
        mle_var = np.zeros(len(self.data))

        if self.config['use_bin_numbers'] == True:
            zbin_col = self.data['bin_number']
        else:
            z_hist = np.histogram(self.data[z_tag], bins=nbins)
            zbin_col = np.digitize(self.data[z_tag], bins=z_hist[1])

        bin_numbers = np.unique(zbin_col)
        for zb in bin_numbers:
            slice = (zbin_col == zb)
            this_est, this_wt = self.optimal_estimator(self.data[slice])
            mle[slice] = this_est
            mle_var[slice] = this_wt

        self.mle = mle
        self.mle_var = mle_var

    def make_treecorr_obj(self):
        '''
        Placeholder for function that can return a treecorr object with k, w if desired!
        '''
        self.tccat = treecorr.Catalog(ra=self.coords.ra.deg,
                            dec=self.coords.dec.deg, ra_units='deg',
                            dec_units='deg', k=redcalc.mle,
                            w=redcalc.mle_var)
