from src.plotter import OverlapPlotter
import os

cat_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_output'
out_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_output'

cat1_basename = 'FgMask_randoms_WISExSCOSmask.fits'
cat2_basename = 'Bg+FgMask_randoms_WISExSCOSmask.fits'
cat3_basename = 'redmagic_hidens_randoms_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'WISExSuperCOSMOS randoms'
cat2_label = 'Masked WISExSuperCOSMOS randoms'
cat3_label = 'redMaGiC hi-dens randoms'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)
cat3 = os.path.join(cat_dir, cat3_basename)

op = OverlapPlotter(cat1_name=cat1)
'''
op2 = OverlapPlotter(cat1_name=cat2)
op3 = OverlapPlotter(cat1_name=cat2, cat2_name=cat3)

outname = os.path.join(out_dir, 'fgmask_scos_randoms_rectilin.png')
op.make_plot(outname=outname, label1=cat1_label, coordframe1='icrs')

outname = os.path.join(out_dir, 'fgmask_scos_hilum_hiz_overlap_aitoff.png')
op.make_plot(outname=outname, projection='aitoff', label1=cat1_label,
             coordframe1='icrs')

outname2 = os.path.join(out_dir, 'double_masked_scos_randoms_rectilin.pdf')
op2.make_plot(outname=outname2, coordframe2='icrs',
              label1=cat2_label)

outname2 = os.path.join(out_dir, 'double_masked_scos_randoms_aitoff.pdf')
op2.make_plot(outname=outname2, projection='aitoff',
             coordframe1='icrs', label1=cat2_label)
'''     
# Make another one
op = OverlapPlotter(cat1_name=cat2, cat2_name=cat3)

outname = os.path.join(out_dir, 'scos_rm_hidens_randoms_overlap_rectilin.pdf')
op.make_plot(outname=outname, label1=cat2_label, coordframe1='icrs', 
             label2=cat3_label, coordframe2='icrs')

outname = os.path.join(out_dir, 'scos_rm_hidens_randoms_overlap_aitoff.pdf')
op.make_plot(outname=outname, label1=cat2_label, coordframe1='icrs',
             label2=cat3_label, coordframe2='icrs', projection='aitoff')
