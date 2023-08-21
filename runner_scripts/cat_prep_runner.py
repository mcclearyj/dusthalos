import os
import pdb
import time
import src.utils as utils
from src.catalog import Catalog
from src.hpmask import HpMask
from src.cat_utils import cat_config_checker

#config = utils.read_yaml('prep_catalog_config.yaml')
config = utils.read_yaml('configs/prep_hidens_randoms_config.yaml')

overwrite = True
vb = True

### TO DO:
###     - Make this a real script
###     - Make error checking more robust?
###     - Specify order of operations in configs?

# maybe a little shady to do this? Ideally, would happen
# somewhere inside Catalog.
config = cat_config_checker(config)

if not os.path.isdir(config['paths']['output_path']):
    os.system(f'mkdir {config["paths"]["output_path"]}')

#run_config = read_yaml(args.config)

# Try loading catalog from the configuration file.
bg_redshift = Catalog(config=config['background_catalog'], vb=vb)

# Load mask from config -- TO DO: add error checking!
bg_redshift.create_mask_from_config(mask_config=config['background_mask'])

# Load other catalog from configuration file
#fg = Catalog(config=config['foreground_catalog'], vb=vb)

# Create its mask
#fg.create_mask_from_config(mask_config=config['foreground_mask'])

# Instatiate background mask
filepath =  os.path.join(config['foreground_mask']['path'], \
                            config['foreground_mask']['filename'])
coordframe = config['foreground_mask']['coordframe']

fg_mask = HpMask(filepath=filepath, coordframe=coordframe)


# Find overlapping masks
bg_redshift.apply_overlapping_masks(mask1=bg_redshift.mask, mask2=fg_mask)


start = time.time()

# Join if requested
if 'match' in config['background_catalog'].keys():
    match_cat = Catalog(config['background_catalog']['match'])
    bg_redshift.match_to_catalog(match_cat, overwrite=overwrite)

end = time.time()
print(f"\n Random matching took {end-start} ms \n")
