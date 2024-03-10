from src.plotter import OverlapPlotter
import os


cat_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs'
out_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'
#cat_dir = '/work/mccleary_group/dusty_halos/prep_cat_output'

cat1_basename = 'prep_cat_output/DoubleMasked_sdss_bg_gals.fits'
cat2_basename = 'prep_cat_output/DoubleMasked_wiseScosPhotoz160708_redshift_cut.fits'

cat1_label = 'SDSS galaxies'
cat2_label = 'WISExSCOS galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(cat1_name=cat1, cat2_name=cat2)

outname = os.path.join(out_dir, 'masked_sdss_wise_overlap_rectilin.png')

op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs', label2=cat2_label,
    coordframe2='galactic', ra_key1='ra', dec_key1='dec', ra_key2='l', dec_key2='b'
)


outname = os.path.join(out_dir, 'masked_sdss_wise_overlap_aitoff.png')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs', label2=cat2_label,
    coordframe2='galactic', ra_key1='ra', dec_key1='dec', ra_key2='l',
    dec_key2='b', projection='mollweide'
)
