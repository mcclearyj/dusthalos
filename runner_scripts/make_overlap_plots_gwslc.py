from src.plotter import OverlapPlotter
import os

cat_dir = '/work/mccleary_group/dusty_halos/catalogs/'
out_dir = '/work/mccleary_group/dusty_halos/catalogs/overlap_plots'

cat1_basename = 'GSWLC-X2_in_SDSS_z_lt_0.18.fits'

cat1_label = 'GWSLC galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)

op = OverlapPlotter(
    cat1_name=cat1, 
    subsample_cat1=False
)

outname = os.path.join(out_dir, 'gwslc_rectilin.pdf')
op.make_plot(
    outname=outname, label1=cat1_label,
    coordframe1='icrs', ra_key1='RA', dec_key1='Dec', 
)

outname = os.path.join(out_dir, 'gwslc_aitoff.pdf')
op.make_plot(
    outname=outname, label1=cat1_label,
    coordframe1='icrs', ra_key1='RA', dec_key1='Dec', 
    projection='aitoff', central_longitude=180
)

# Save memory!
del(op)



###
### Also do it for randoms!
###

cat1_basename = 'gwslc_randoms.fits'

cat1_label = 'GWSLC randoms'

cat1 = os.path.join(cat_dir, cat1_basename)

# Randoms are large, so just plot a subsample
op = OverlapPlotter(
    cat1_name=cat1,
    subsample_cat1=True, subsample_size=3e5
)

outname = os.path.join(out_dir, 'rand_gwslc_rectilin.pdf')
op.make_plot(
    outname=outname, label1=cat1_label,
    coordframe1='icrs', ra_key1='ra', dec_key1='dec', 
)

outname = os.path.join(out_dir, 'rand_gwslc_aitoff.pdf')
op.make_plot(
    outname=outname, label1=cat1_label,
    coordframe1='icrs', ra_key1='ra', dec_key1='dec', 
    central_longitude=180, projection='aitoff'
)
