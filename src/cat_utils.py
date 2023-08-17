import os
from astropy.table import Table, join
import numpy as np
from . import utils

def cat_checker(catalog):
    '''
    Minimal helper function to ensure required attributes are present. We could
    think about making this a cat *loader* too, not just a cat checker
    '''
    # A dict is fine too
    if type(catalog) == dict:
        catalog = utils.AttrDict(catalog)

    required_attributes = ['filepath', 'tabname', 'data']

    for key in required_attributes:
        if key not in catalog.__dict__.keys():
            raise KeyError(f'Catalog object is missing attribute {key}')

def cat_config_checker(config):
    '''
    Utility method to check
    '''

    cf_keys = config.keys()
    # Required catalog keys
    required_cat_keys = ['filename', 'ra_tag', 'dec_tag',
                        'coordframe', 'coord_units']
    required_match_keys = ['filename', 'match_type', 'tabname']
    # Required mask keys
    required_mask_keys = ['filename', 'coordframe']

    # Go through and ensure that minimal required keywords are present
    if 'paths' not in cf_keys:
        raise KeyError('config file is missing "paths" setting')
    if 'catalog1' not in cf_keys:
        raise KeyError('config file is missing "catalog1"')

    # Make sure there is a path defined
    for sub_cfg in cf_keys:
        if ('path' not in config[sub_cfg].keys()) or \
            (config[sub_cfg]['path'] is None):
                config[sub_cfg]['path'] = config['paths']['catalog_path']
        if ('output_path' not in config[sub_cfg].keys()) or \
            (config[sub_cfg]['output_path'] is None):
                config[sub_cfg]['output_path'] = config['paths']['output_path']

    # For catalog1, catalog2: make sure the coord tags are there
    catalogs = [config['catalog1']]
    masks = []

    if 'catalog2' in cf_keys:
        catalogs.append(config['catalog2'])

    if 'catalog1_mask' in cf_keys:
        masks.append(config['catalog1_mask'])

    if 'catalog2_mask' in cf_keys:
        masks.append(config['catalog2_mask'])

    for cat in catalogs:
        cat_cfg = cat.keys()
        for key in required_cat_keys:
            if key not in cat_cfg:
                raise KeyError(f'catalog config is missing required key: {key}')
            if key == 'coord_units':
                if cat[key] not in ['deg', 'rad']:
                    print(f'The "coord_units" parameter supplied was {key} ' +
                        'but must be either "deg" or "rad"')
        if ('match' in cat.keys()) & (cat['match'] != None):
            for mkey in required_match_keys:
                if mkey not in cat['match']:
                    raise KeyError(f'match config missing key: {mkey}')

    for mask in masks:
        mask_cfg = cat.keys()
        for key in required_mask_keys:
            if key not in mask_cfg:
                raise KeyError(f'catalog config is missing required key: {key}')

    return config
