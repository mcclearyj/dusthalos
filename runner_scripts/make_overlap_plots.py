from src.plotter import OverlapPlotter
import os


cat_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs'
out_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'

cat1_basename = 'prep_cat_output/DoubleMasked_sdss_bg_randoms.fits'
cat2_basename = 'prep_cat_output/DoubleMasked_sdss_bg_gals.fits'
cat3_basename = 'gaia_desdr2_overlap_full.fits'

cat1_label = 'Randoms'
cat2_label = 'Galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)
cat3 = os.path.join(cat_dir, cat3_basename)

op = OverlapPlotter(cat1_name=cat2, cat2_name=cat1)
op2 = OverlapPlotter(cat1_name=cat3, cat2_name=cat2)

outname = os.path.join(out_dir, 'masked_sdss_bg_overlap_rectilin.pdf')
op.make_plot(outname=outname, label1=cat2_label, coordframe1='icrs',
             label2=cat1_label, coordframe2='icrs',
             ra_key1='ra_rand', dec_key1='dec_rand'
             )

outname = os.path.join(out_dir, 'masked_sdss_bg_overlap_aitoff.pdf')
op.make_plot(outname=outname, projection='aitoff', label1=cat2_label,
             coordframe1='icrs', label2=cat1_label, coordframe2='icrs',
             ra_key1='ra_rand', dec_key1='dec_rand'
             )
'''
outname2 = os.path.join(out_dir, 'gaia_redmagic_overlap_rectilin.pdf')
op2.make_plot(outname=outname, projection='aitoff', label1=cat1_label,
             coordframe1='icrs', label2=cat3_label, coordframe2='icrs',
             ra_key2='ra_redmagic_hiz', dec_key2='dec_redmagic_hiz'
             )

outname2 = os.path.join(out_dir, 'gaia_redmagic_overlap_aitoff.pdf')
op2.make_plot(outname=outname, projection='aitoff', label1=cat1_label,
             coordframe1='icrs', label2=cat3_label, coordframe2='icrs',
             ra_key2='ra_redmagic_hiz', dec_key2='dec_redmagic_hiz'
             )
'''
