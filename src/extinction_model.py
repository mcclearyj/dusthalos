import numpy as np
import pdb
import extinction
import numpy as np
from dust_extinction.parameter_averages import CCM89, G23
import astropy.units as u

class ExtinctionModel:
    """
    Wrapper for the extinction python library
    """
    def __init__(self, dmconfig):
        self.dmdp = []
        self.dmconfig = dmconfig
        self.allowed_models = [
            'ccm89', 'odonnell94', 'calzetti00',
            'fitzpatrick99', 'g23', 'g16']

    def _dust_config_checker(self):
        """
        Make sure minimal configuration parameters are present and also ensure
        that given model is allowed
        """
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
        """
        Get model
        """
        r_v = self.dmconfig['R']
        model = self.dmconfig['model']
        kwargs = {
            'wave': np.array(self.dmconfig['wavelengths'])*u.AA,
            'a_v': self.dmconfig['A_V']
        }
        if model == 'ccm89':
            #self.dmdp = extinction.ccm89(**kwargs, r_v=r_v)
            ccm89 = CCM89(Rv=r_v)
            self.dmdp = ccm89.evaluate(kwargs['wave'], Rv=r_v)
        elif model == 'odonnell94':
            self.dmdp = extinction.odonnell94(**kwargs, r_v=r_v)
        elif model == 'fitz99':
            self.dmdp = extinction.fitzpatrick99(**kwargs, r_v=r_v)
        elif model == 'calzetti00':
            self.dmdp = extinction.calzetti00(**kwargs, r_v=r_v)
            # array([1.65709158, 1.22587063, 0.84134317, 0.59816654, 0.40424546])
        elif model == 'g23':
            gordon23 = G23(Rv=r_v)
            self.dmdp = gordon23.evaluate(kwargs['wave'], Rv=r_v)
            # array([1.55892008, 1.20243405, 0.85918173, 0.6445867, 0.4829834])

        # Comfort display
        print(f"\nUsing extinction model {model}: {self.dmdp}\n")

    def get_dust_model(self):
        """
        Get a dust model based on configuration file parameters
        """
        # Check configuration file
        self._dust_config_checker()

        # Get model
        self._get_dust_model()
