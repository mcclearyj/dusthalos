from src.plotter import OverlapPlotter
import os


cat_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs'
out_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'

cat1_basename = 'prep_cat_output/DoubleMasked_gaia_desdr2_overlap_full.fits'
cat2_basename = 'prep_cat_output/redmagic_hiz_y3_GOLD_JOINED_catalog.fits'
cat3_basename = 'gaia_desdr2_overlap_full.fits'

cat1_label = 'Masked Gaia stars'
cat2_label = 'redMaGiC hi-lum hi-z'
cat3_label = 'Gaia stars'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)
cat3 = os.path.join(cat_dir, cat3_basename)

op = OverlapPlotter(cat1_name=cat2, cat2_name=cat1)
op2 = OverlapPlotter(cat1_name=cat3, cat2_name=cat2)

outname = os.path.join(out_dir, 'masked_gaia_redmagic_overlap_rectilin.pdf')
op.make_plot(outname=outname, label1=cat2_label, coordframe1='icrs',
             label2=cat1_label, coordframe2='icrs',
             ra_tag1='ra_redmagic_hiz', dec_tag1='dec_redmagic_hiz'
             )

outname = os.path.join(out_dir, 'masked_gaia_redmagic_overlap_aitoff.pdf')
op.make_plot(outname=outname, projection='aitoff', label1=cat2_label,
             coordframe1='icrs', label2=cat1_label, coordframe2='icrs',
             ra_tag1='ra_redmagic_hiz', dec_tag1='dec_redmagic_hiz'
             )
'''
outname2 = os.path.join(out_dir, 'gaia_redmagic_overlap_rectilin.pdf')
op2.make_plot(outname=outname, projection='aitoff', label1=cat1_label,
             coordframe1='icrs', label2=cat3_label, coordframe2='icrs',
             ra_tag2='ra_redmagic_hiz', dec_tag2='dec_redmagic_hiz'
             )

outname2 = os.path.join(out_dir, 'gaia_redmagic_overlap_aitoff.pdf')
op2.make_plot(outname=outname, projection='aitoff', label1=cat1_label,
             coordframe1='icrs', label2=cat3_label, coordframe2='icrs',
             ra_tag2='ra_redmagic_hiz', dec_tag2='dec_redmagic_hiz'
             )
'''
