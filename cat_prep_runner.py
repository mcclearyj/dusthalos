import os
import pdb

import src.utils as utils
from src.catalog import Catalog

# First, match the background photometry and redshift catalogs
# In very short order, this needs to be a config of some sort...

catalog_path = '/Users/j.mccleary/Research/dusty_halos/catalogs'
redshift_name = 'y3a2_gold2.2.1_redmagic_highlum_highz.fits'
bgmask_name = 'y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits'
phot_cat_name = '14741.parquet'
fg_cat_mask = 'WISExSCOSmask.fits'


bg_redshift = Catalog(catfile='y3a2_gold2.2.1_redmagic_highlum_highz.fits',
                        datadir=catalog_path, tabname='z_rm')

bg_redshift.join_cats_kw(cat2=os.path.join(catalog_path, '14742.parquet'),
                            name2='z_rm')
