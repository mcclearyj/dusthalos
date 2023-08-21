import os
import pdb
import time
import src.utils as utils
from src.hpmask import HpMask
from src.fg_randoms import FgRandoms
from src.cat_utils import cat_config_checker

config = utils.read_yaml('configs/prep_randoms.yaml')
overwrite = True
vb = True

# maybe a little shady to do this?
config = cat_config_checker(config)

# Instantiate FgRandoms.
fgr = FgRandoms(config=config['foreground_mask'])

# Make randoms on the sphere
fgr.make_fg_randoms(overwrite=overwrite)

# Instatiate background mask
filepath =  os.path.join(config['background_mask']['path'], \
                            config['background_mask']['filename'])
coordframe = config['background_mask']['coordframe']

bg_mask = HpMask(filepath=filepath, coordframe=coordframe)

# Filter coordinates according to background mask specified in config
fgr.do_overlapping_mask(bg_mask, overwrite=overwrite)
