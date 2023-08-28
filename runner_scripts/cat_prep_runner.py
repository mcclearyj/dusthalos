import os
import time
import argparse
import src.utils as utils
from src.catalog import Catalog
from src.cat_utils import cat_config_checker

def main(args):
    config_path = args.config_path
    config = utils.read_yaml(config_path)

    overwrite = True
    vb = True

    # Checking configuration consistency
    config = cat_config_checker(config)

    # Create output directory if it doesn't exist
    if not os.path.isdir(config['paths']['output_path']):
        os.makedirs(config['paths']['output_path'])

    # Load foreground catalog from configuration file
    fg = Catalog(config=config['foreground_catalog'], vb=vb)

    # Create its mask
    fg.create_mask_from_config(mask_config=config['foreground_mask'])

    # Apply foreground mask
    fg.appy_mask()
    # Load background catalog from the configuration file
    bg_redshift = Catalog(config=config['background_catalog'], vb=vb)

    # Load mask from config
    bg_redshift.create_mask_from_config(mask_config=config['background_mask'])

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
    parser.add_argument('--config_path', type=str, help='Path to the configuration file.', required=True)
    args = parser.parse_args()
    main(args)
