import os
import pdb
import time
import src.utils as utils
from src.catalog import Catalog
from src.cat_utils import cat_config_checker

config = utils.read_yaml('prep_catalog_config.yaml')
overwrite = True
vb = True

# maybe a little shady to do this?
config = cat_config_checker(config)

if not os.path.isdir(config['paths']['output_path']):
    os.system(f'mkdir {config["paths"]["output_path"]}')

#run_config = read_yaml(args.config)

# Try loading catalog from the configuration file.
bg_redshift = Catalog(config=config['catalog1'], vb=vb)

# Load mask from config -- TO DO: add error checking!
bg_redshift.create_mask_from_config(mask_config=config['catalog1_mask'])

# Load other catalog from configuration file
fg = Catalog(config=config['catalog2'], vb=vb)

# Create its mask
fg.create_mask_from_config(mask_config=config['catalog2_mask'])

# Find overlapping masks
bg_redshift.apply_overlapping_masks(mask1=bg_redshift.mask, mask2=fg.mask)

# Join if requested
if 'match' in config['catalog1'].keys():
    match_cat = Catalog(config['catalog1']['match'])
    bg_redshift.match_to_catalog(match_cat, overwrite=overwrite)
