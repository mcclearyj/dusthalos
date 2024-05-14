from src.plotter import OverlapPlotter
import os

cat_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_gwslc'
out_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_gwslc/overlap_plots'

cat1_basename = 'sf/DoubleMasked_sf_gswlc_galaxies.fits'
cat2_basename = 'DoubleMasked_sdss_bg_photoz.fits'

cat1_label = 'GWSLC star-forming galaxies'
cat2_label = 'SDSS galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(
    cat1_name=cat1, cat2_name=cat2, 
    subsample_cat1=False, subsample_cat2=True, subsample_size=1e6
)


outname = os.path.join(out_dir, 'masked_gwslc_sf_sdss_overlap_rectilin.png')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs', label2=cat2_label,
    coordframe2='icrs', ra_key1='RA', dec_key1='Dec', 
    ra_key2='ra', dec_key2='dec', central_longitude=0
)

outname = os.path.join(out_dir, 'masked_gwslc_sf_sdss_overlap_aitoff.png')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs', label2=cat2_label,
    coordframe2='icrs', ra_key1='RA', dec_key1='Dec', 
    ra_key2='ra', dec_key2='dec', projection='aitoff', central_longitude=180
)

# Save memory!
del(op)



###
### Also do it for randoms!
###

cat1_basename = 'sf/DoubleMasked_sf_gwslc_randoms.fits'
cat2_basename = 'rand_sdss_bg_JOINED_catalog.fits'

cat1_label = 'GWSLC star-forming randoms'
cat2_label = 'SDSS randoms'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

# Randoms are large, so just plot a subsample
op = OverlapPlotter(
    cat1_name=cat1, cat2_name=cat2, 
    subsample_cat1=True, subsample_cat2=True, subsample_size=1e6
)

outname = os.path.join(out_dir, 'rand_gwslc_sf_sdss_overlap_rectilin.png')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs', label2=cat2_label,
    coordframe2='icrs', ra_key1='ra', dec_key1='dec', 
    ra_key2='ra_rand', dec_key2='dec_rand', central_longitude=0
)

outname = os.path.join(out_dir, 'rand_gwslc_sf_sdss_overlap_aitoff.png')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs', label2=cat2_label,
    coordframe2='icrs', ra_key1='ra', dec_key1='dec',
    ra_key2='ra_rand', dec_key2='dec_rand', central_longitude=180,
    projection='aitoff',
)
