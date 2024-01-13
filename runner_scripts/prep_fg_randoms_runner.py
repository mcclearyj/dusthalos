import os
import pdb
import time
import argparse
import src.utils as utils
from src.hpmask import HpMask
from src.fg_randoms import FgRandoms
from src.cat_utils import all_config_checker


def main(args):
    
    config_file = args.config_file
    config = utils.read_yaml(config_file)
    
    overwrite = True
    vb = True

    # maybe a little shady to do this?
    config = all_config_checker(config)

    # Instantiate FgRandoms.
    fgr = FgRandoms(config=config['foreground_mask'])
    
    # Make randoms on the sphere
    fgr.make_fg_randoms(nrand=3e7)
    
    # Instatiate background mask
    filepath =  os.path.join(config['background_mask']['path'], \
                                config['background_mask']['filename'])
    coordframe = config['background_mask']['coordframe']
    
    bg_mask = HpMask(filepath=filepath, coordframe=coordframe)
    
    # Filter coordinates according to background mask specified in config
    fgr.do_overlapping_mask(bg_mask, overwrite=overwrite)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
             description='Runner script for Catalog operations.')
    parser.add_argument('--config_file', type=str, 
                        help='Configuration file.', required=True)
    args = parser.parse_args()
    main(args)
    