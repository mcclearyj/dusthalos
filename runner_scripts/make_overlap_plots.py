from src.plotter import OverlapPlotter
import os


#cat_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs'
#out_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
cat_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_output'
out_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_output/overlap_plots'

cat1_basename = 'DoubleMasked_wiseScosPhotoz160708_redshift_cut.fits'
cat2_basename = 'redmagic_hidens_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'WISExSCOS galaxies'
cat2_label = 'redMaGiC hi-dens galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(cat1_name=cat1, cat2_name=cat2)

outname = os.path.join(out_dir, 'masked_scos_rmz_hidens_overlap_rectilin.png')

op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='galactic', label2=cat2_label,
    coordframe2='icrs', ra_key1='l', dec_key1='b', 
    ra_key2='ra_redmagic_hidens', dec_key2='dec_redmagic_hidens'
)


outname = os.path.join(out_dir, 'masked_scos_rmz_hidens_overlap_aitoff.png')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='galactic', label2=cat2_label,
    coordframe2='icrs', ra_key1='l', dec_key1='b', 
    ra_key2='ra_redmagic_hidens', dec_key2='dec_redmagic_hidens', projection='aitoff'
)

###
### Also do it for randoms!
###

cat2_basename = 'DoubleMasked_WISExSCOS_randoms3.fits'
cat1_basename = 'redmagic_hidens_randoms_y3_GOLD_JOINED_catalog.fits'

cat2_label = 'WISExSCOS randoms'
cat1_label = 'redMaGiC hi-dens randoms'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(cat1_name=cat1, cat2_name=cat2)

outname = os.path.join(out_dir, 'rand_scos_rmz_hidens_overlap_rectilin_v2.png')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs', label2=cat2_label,
    coordframe2='icrs', ra_key2='ra', dec_key2='dec', 
    ra_key1='ra_redmagic_hidens_randoms', dec_key1='dec_redmagic_hidens_randoms'
)

outname = os.path.join(out_dir, 'rand_scos_rectilin_v2.png')
op.make_plot(
    outname=outname, label1=cat2_label, coordframe1='icrs', ra_key1='ra', dec_key1='dec',
)

