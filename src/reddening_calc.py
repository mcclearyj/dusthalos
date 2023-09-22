import numpy as np
import pdb
from src.extinction_model import ExtinctionModel

class ReddeningCalculator(ExtinctionModel):

    def __init__(self, data, redcalc_config=None):
        '''
        Do reddening calculation
        '''
        self.data = data # Data with magnitudes, redshifts
        self.mle = [] # Maximum likelihood estimator for anomalous reddening
        self.mle_var = [] # Variance of said ML
        self.nbins = None
        self.z_tag = None
        self.use_bin_numbers = False

        # dust config should be in here now too
        self.load_config(redcalc_config)

        super().__init__(dmconfig = redcalc_config['dust_params'])

    def load_config(self, redcalc_config):
        '''
        redshift config should have all the parameters keywords
        '''
        # Resonable redshift defaults, not sure about dust params
        #default_dust_pars = {'model': 'calzetti00', 'R': 3.1, 'A_V': 1,
        #                    'wavelengths': [4808.49, 6417.65, 7814.58, 9168.85]}
        # No, it should fail loudly.

        base_config = {'nbins': 7,
                          'z_tag': 'z',
                          'use_bin_numbers': False,
                          'dust_params': None,
                          'bandpasses': None,
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

        print(f'\n ReddeningCalculator configuration values:')
        print(f'\t {base_config}')

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
            band_bool = (band_mag > -9999) & \
                        (band_mag != np.nan) & \
                        (band_mag < 30)
            wg *= band_bool

        # percent of galaxies that failed
        pfail = 100-(np.count_nonzero(wg) / catlen * 100)

        # How many failed?
        print(f'ReddeningCalc: {np.count_nonzero(wg)}/{catlen} galaxies',
              f'({100-pfail:.1f}%) have good photometry')
        print(f'ReddeningCalc: removing {pfail:.1f}% of galaxies from data')
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
        dmdp = self.dmdp

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

    def run(self):
        '''
        zmin: minimum redshift of galaxies to consider (default None)
        nbins: number of bins for color estimation
        '''

        print(f'Removing catalog entries with NaN & sentinel values')
        self.preprocess_catalog()

        print(f'Make dust extinction model')
        self.get_dust_model()

        print(f'Beginning background catalog reddening calculation...')

        mle = np.zeros(len(self.data))
        mle_var = np.zeros(len(self.data))

        if self.redcalc_config['use_bin_numbers'] == True:
            zbin_col = self.data['bin_number']
        else:
            # Make a "bin number" on the fly!
            z_hist = np.histogram(self.data[self.redcalc_config['z_tag']],
                                    bins=self.redcalc_config['nbins'])
            zbin_col = np.digitize(self.data[self.redcalc_config['z_tag']],
                                    bins=z_hist[1])

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
        self.treecorr_cat = treecorr.Catalog(ra=self.coords.ra.deg,
                            dec=self.coords.dec.deg, ra_units='deg',
                            dec_units='deg', k=redcalc.mle,
                            w=redcalc.mle_var)
