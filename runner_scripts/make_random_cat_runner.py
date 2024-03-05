import os
import pdb
import time
import argparse
import src.utils as utils
from src.hpmask import HpMask
from src.random_cat import RandomCat
from src.cat_utils import all_config_checker


def main(args):
    #setting overwrite=True; may make configurable later
    overwrite = True

    # Load config file
    config_file = args.config
    run_config = utils.read_yaml(config_file)

    # This step populates default paths too -- may be a bit shady to do this
    run_config = all_config_checker(run_config)

    # Instantiate FgRandoms.
    fgr = RandomCat(
        mask_config=run_config['mask1'], nrand=run_config.get('nrand'),
        output_name=run_config.get('output_name')
    )

    # Make randoms on the sphere
    fgr.make_random_cat()

    # Filter random points against mask2 if supplied in config
    if run_config.get('mask2') != None:
        filepath = os.path.join(
            run_config['mask2']['path'], run_config['mask2']['filename']
        )

        mask2 = HpMask(
            filepath=filepath, coordframe=run_config['mask2']['coordframe']
        )

        fgr.do_overlapping_mask(
            mask2, catalog_for_comparison=run_config.get('comparison_cat'),
            overwrite=overwrite
        )

if __name__ == '__main__':

    # Load arguments
    parser = argparse.ArgumentParser(
        description='Runner script for Catalog operations.'
    )
    parser.add_argument(
        '-config', '-c', type=str, help='Configuration file', required=True
    )
    args = parser.parse_args()

    # Go!
    main(args)
