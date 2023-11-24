from src.plotter import OverlapPlotter
import os

cat_dir = '/Users/j.mccleary/Research/dusty_halos/dcatalogs/prep_cat_output'
out_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'

cat1_basename = 'FgMask_randoms_WISExSCOSmask.fits'
cat2_basename = 'Bg+FgMask_randoms_WISExSCOSmask.fits'
cat3_basename = 'redmagic_hiz_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'WISExSuperCOSMOS randoms'
cat2_label = 'Masked WISExSuperCOSMOS randoms'
cat3_label = 'redMaGiC hi-lum hi-z'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)
cat3 = os.path.join(cat_dir, cat3_basename)

op = OverlapPlotter(cat1_name=cat1)
op2 = OverlapPlotter(cat1_name=cat2)

outname = os.path.join(out_dir, 'full_scos_randoms_rectilin.png')
op.make_plot(outname=outname, label1=cat1_label, coordframe1='icrs')

outname = os.path.join(out_dir, 'full_scos_hilum_hiz_overlap_aitoff.png')
op.make_plot(outname=outname, projection='aitoff', label1=cat1_label,
             coordframe1='icrs')


outname2 = os.path.join(out_dir, 'masked_scos_randoms_rectilin.pdf')
op2.make_plot(outname=outname2, coordframe2='icrs',
              label1=cat2_label)

outname2 = os.path.join(out_dir, 'masked_scos_randoms_aitoff.pdf')
op2.make_plot(outname=outname2, projection='aitoff',
             coordframe1='icrs', label1=cat2_label)
