import numpy as np
import pdb
import extinction
import numpy as np

class ExtinctionModel:
    '''
    Wrapper for the extinction python library
    '''
    def __init__(self, dmconfig):
        self.dmdp = []
        self.dmconfig = dmconfig
        self.allowed_models = ['ccm89', 'odonnell94', 'calzetti00',
                               'fitzpatrick99', 'fm07']

    def _dust_config_checker(self):
        '''
        Make sure minimal configuration parameters are present and also ensure
        that given model is allowed
        '''
        required_keys = ['model', 'R', 'A_V', 'wavelengths']
        allowed = self.allowed_models
        config = self.dmconfig

        cfg_keys = config.keys()

        for key in required_keys:
            if key not in cfg_keys:
                msg = f'\nExtinction: config missing required key {key}\n'
                raise ValueError(msg)

        if config['model'] not in allowed:
                msg = f'\nExtinction: {config["model"]} not in {allowed}'
                raise NameError(msg)

    def _get_dust_model(self):
        '''
        Get model
        '''
        r_v = self.dmconfig['R']
        model = self.dmconfig['model']
        kwargs = {
            'wave': np.array(self.dmconfig['wavelengths']),
            'a_v': self.dmconfig['A_V']
        }
        if model == 'ccm89':
            self.dmdp = extinction.ccm89(**kwargs, r_v=r_v)
        elif model == 'odonnell94':
            self.dmdp = extinction.odonnell94(**kwargs, r_v=r_v)
        elif model == 'fitz99':
            self.dmdp = extinction.fitz99(**kwargs, r_v=r_v)
        elif model == 'calzetti00':
            self.dmdp = extinction.calzetti00(**kwargs, r_v=r_v)
        else:
            model = 'fm07'
            self.dmdp = extinction.fm07(**kwargs)

        print(f"Using extinction model {model}") # Comfort display

    def get_dust_model(self):
        '''
        Get a dust model based on configuration file parameters
        '''
        # Check configuration file
        self._dust_config_checker()

        # Get model
        self._get_dust_model()
