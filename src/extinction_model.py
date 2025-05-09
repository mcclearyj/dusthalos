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
        # Note: dir() returns everything, not just classes, but OK for now
        self.parameter_averages_models = dir(parameter_averages)
        self.averages_models = dir(averages)

    def _get_allowed_models(self):
        """ 
        Placeholder for function that would do model-checking in a more robust way 
        """
        parameter_averages_models = [
            i[0] for i in inspect.getmembers(parameter_averages, inspect.isclass)
        ]
        averages_models = [
            i[0] for i in inspect.getmembers(averages, inspect.isclass)
        ]

        allowed_models = []
        allowed_models.extend(parameter_averages_models)
        allowed_models.extend(averages_models)
        
        return allowed_models 
        

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

        # This try/except is a bit redundant, since config checker catches this first.
        if (model in self.parameter_averages_models) | (model in self.averages_models):
            if model in self.parameter_averages_models:
                dust_method = getattr(parameter_averages, model)()
            else:
                dust_method = getattr(averages, model)()
            if 'Rv' in inspect.signature(dust_method.evaluate).parameters.keys():
                self.dmdp = dust_method.evaluate(kwargs['wave'], Rv=r_v)
            else:
                self.dmdp = dust_method.evaluate(kwargs['wave'])
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
