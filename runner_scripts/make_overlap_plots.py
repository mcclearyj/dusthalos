from src.plotter import OverlapPlotter
import os

cat_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
outdir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
cat1_basename = 'DoubleMasked_Masked_wiseScosPhotoz160708_redshift_cut.fits'
cat2_basename = 'redmagic_hiz_y3_GOLD_JOINED_catalog.fits'
cat3_basename = 'redmagic_hidens_y3_GOLD_JOINED_catalog.fits'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)
cat3 = os.path.join(cat_dir, cat3_basename)

op = OverlapPlotter(cat1_name=cat1, cat2_name=cat2)
op2 = OverlapPlotter(cat1_name=cat1, cat2_name=cat3)

outname = os.path.join(outdir, 'masked_hilum_hiz_rectilin.pdf')
op.make_plot(outname=outname, ra_tag2='ra_redmagic_hiz', dec_tag2='dec_redmagic_hiz', coordframe2='icrs')

outname = os.path.join(outdir, 'masked_hilum_hiz_mollweide.pdf')
op.make_plot(outname=outname, projection='mollweide', ra_tag2='ra_redmagic_hiz',
             dec_tag2='dec_redmagic_hiz', coordframe2='icrs')


outname2 = os.path.join(outdir, 'masked_hidens_rectilin.pdf')
op2.make_plot(outname=outname2, ra_tag2='ra_redmagic_hidens',
                     dec_tag2='dec_redmagic_hidens', coordframe2='icrs')

outname2 = os.path.join(outdir, 'masked_hidens_mollweide.pdf')
op2.make_plot(outname=outname2, projection='mollweide', ra_tag2='ra_redmagic_hidens',
                     dec_tag2='dec_redmagic_hidens', coordframe2='icrs')
