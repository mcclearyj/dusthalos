import numpy as np
import pdb
import inspect
from dust_extinction import averages, parameter_averages
import astropy.units as u

class ExtinctionModel:
    """
    Wrapper for the extinction python library
    """
    def __init__(self, dmconfig):
        self.dmdp = []
        self.dmconfig = dmconfig
        # Adds a bunch of other stuff, not just models, but OK for now.
        self.allowed_deprecated_models = dir(extinction)
        self.averages_models = dir(averages)
        self.parameter_averages_models = dir(parameter_averages)


    def _dust_config_checker(self):
        """
        Make sure minimal configuration parameters are present and also ensure
        that given model is allowed
        """
        required_keys = ['model', 'R', 'A_V', 'wavelengths']
        config = self.dmconfig
        cfg_keys = config.keys()

        for key in required_keys:
            if key not in cfg_keys:
                msg = f'\nExtinction: config missing required key {key}\n'
                raise ValueError(msg)
                
        model = config['model']
        if not (hasattr(parameter_averages, model) or hasattr(averages, model)):
            msg = f'\nExtinction: {model} not found in parameter_averages or extinction\n'
            raise NameError(msg)

    """
    def _dust_config_checker(self):
        #Make sure minimal configuration parameters are present and also ensure
        #that given model is allowed
        required_keys = ['model', 'R', 'A_V', 'wavelengths']
        allowed = self.allowed_deprecated_models
        allowed.append(self.averages_models + self.parameter_averages_models)
        config = self.dmconfig

        cfg_keys = config.keys()

        for key in required_keys:
            if key not in cfg_keys:
                msg = f'\nExtinction: config missing required key {key}\n'
                raise ValueError(msg)

	pdb.set_trace()

        if config['model'] not in allowed:
                msg = f'\nExtinction: {config["model"]} not in {allowed}'
                raise NameError(msg)
    """
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

        if (model in self.parameter_averages_models) | (model in self.averages_models):
            if model in self.parameter_averages_models:
                dust_method = getattr(parameter_averages, model)()
            else:
                dust_method = getattr(averages, model)()
            if 'Rv' in inspect.signature(dust_method.evaluate).parameters.items():
                self.dmdp = dust_method.evaluate(kwargs['wave'], Rv=r_v)
            else:
                self.dmdp = dust_method.evaluate(kwargs['wave'])
        elif model in self.allowed_deprecated_models:
            dust_method = getattr(extinction, model)
            self.dmdp = dust_method(**kwargs, r_v=r_v)
        else:
            msg = "Problem selecting dust extinction model, exiting"
            raise ValueError(msg)

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
