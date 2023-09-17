'''
This script runs through the catalog configuration file and
creates masked (and matched, if set to True) catalogs ready for dust
correlation function calculations.

TO DO
    - Allow creation of just foreground or just background catalogs.
      Perhaps a run config could work?
    - Make this a proper Runner class -> could be helpful for allowing
      individual catalog calculations.
'''

import os
import time
import argparse
import src.utils as utils
from src.catalog import Catalog
from src.cat_utils import all_config_checker

def main(args):
    config_file = args.config_file
    config = utils.read_yaml(config_file)

    overwrite = True
    vb = True

    # Checking configuration consistency
    config = all_config_checker(config)

    # Create output directory if it doesn't exist
    if not os.path.isdir(config['paths']['output_path']):
        os.makedirs(config['paths']['output_path'])

    # Load foreground catalog from configuration file
    fg = Catalog(config=config['foreground_catalog'], vb=vb)

    # Create foreground catalog mask from configuration file
    fg.create_mask_from_config(mask_config=config['foreground_mask'])

    # Load background catalog from the configuration file
    bg_redshift = Catalog(config=config['background_catalog'], vb=vb)

    # Create background catalog configuration file
    bg_redshift.create_mask_from_config(mask_config=config['background_mask'])

    # Apply overlapping masks
    #fg.apply_overlapping_masks(mask1=bg_redshift.mask, mask2=fg.mask)

    # Find overlapping masks
    bg_redshift.apply_overlapping_masks(mask1=bg_redshift.mask, mask2=fg.mask)

    start = time.time()

    # Join if requested
    if 'match' in config['background_catalog'].keys():
        match_cat = Catalog(config['background_catalog']['match'], memmap=True)
        bg_redshift.match_to_catalog(match_cat, overwrite=overwrite)

    end = time.time()
    print(f"\n Random matching took {((end-start)/60.):.1f} mins \n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runner script for Catalog operations.')
    parser.add_argument('--config_file', type=str, help='Configuration file.', required=True)
    args = parser.parse_args()
    main(args)
