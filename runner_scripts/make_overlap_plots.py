from src.plotter import OverlapPlotter
import os

cat_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
outdir = '/Users/j.mccleary/Research/dusty_halos/catalogs/overlap_plots'
cat1_basename = 'masked_y3a2_gold2.2.1_redmagic_highlum_highz.fits'
#cat2_basename = 'y3a2_gold2.2.1_redmagic_highlum_highz.fits'

cat1 = os.path.join(cat_dir, cat1_basename)
#cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(cat1_name=cat1) #, cat2_name=cat2)

outname = os.path.join(outdir, 'masked_hilum_hiz_rectilin.pdf')
op.make_plot(outname=outname)

outname = os.path.join(outdir, 'masked_hilum_hiz_mollweide.pdf')
op.make_plot(outname=outname, projection='mollweide')
