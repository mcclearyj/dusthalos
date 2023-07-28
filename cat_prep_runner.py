import .utils
from .catalog_prep import BkgCat, FgRandoms
import os

# First, match the background photometry and redshift catalogs
# In very short order, this needs to be a config of some sort...

catalog_path = '/Users/j.mccleary/Research/dusty_halos/catalogs'
redshift_name = 'y3a2_gold2.2.1_redmagic_highlum_highz.fits'
bgmask_name = 'y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits'
phot_cat_name = '14741.parquet'
fg_cat_mask = 'WISExSCOSmask.fits'

'''
def match_coords(cat1, cat2, ratag1=None, dectag1=None,
                ratag2=None, dectag2=None, radius=0.5):
'''

# Doing a join would be more solid
bg_phot = Table.read(os.path.join(catalog_path, phot_cat_name))
bg_redshift = Table.read(os.path.join(catalog_path, redshift_name))
