from src.plotter import OverlapPlotter
import os

cat_dir = '/n23data1/mccleary/dustyhalos/dusthalos'
out_dir = '/n23data1/mccleary/dustyhalos/dusthalos/overlap_plots'

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
cat1_basename = 'GSWLC-X2_in_SDSS_randoms.fits'
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
### Save memory!
del(op)


###
### Now do cross-overlap for SDSS and GSWLC
###
cat_dir = '/n23data1/mccleary/dustyhalos/dusthalos/prep_cat_gswlc'
out_dir = '/n23data1/mccleary/dustyhalos/dusthalos/prep_cat_gswlc/overlap_plots'

cat1_basename = 'DoubleMasked_GSWLC-X2_in_SDSS_z_lt_0.18.fits'
cat2_basename = 'DoubleMasked_sdss_bg_photoz2.fits'
cat1_label = 'GWSLC galaxies'
cat2_label = 'SDSS galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(
    cat1_name=cat1, 
    cat2_name=cat2,
    subsample_cat1=False,
    subsample_cat2=True, 
    subsample_size=3e5
)

outname = os.path.join(out_dir, 'gwslc_x_sdss_rectilin.png')
op.make_plot(
    outname=outname, label1=cat1_label, label2=cat2_label,
    coordframe1='icrs', ra_key1='RA', dec_key1='Dec', 
    coordframe2='icrs', ra_key2='ra', dec_key2='dec',
)

outname = os.path.join(out_dir, 'gwslc_x_sdss_aitoff.png')
op.make_plot(
    outname=outname, label1=cat1_label, label2=cat2_label,
    coordframe1='icrs', ra_key1='RA', dec_key1='Dec', 
    coordframe2='icrs', ra_key2='ra', dec_key2='dec',
    projection='aitoff', central_longitude=180
)

# Save memory!
del(op)

###
### Finally do cross-overlap for SDSS and GSWLC randoms
###
cat_dir = '/n23data1/mccleary/dustyhalos/dusthalos/prep_cat_gswlc'
out_dir = '/n23data1/mccleary/dustyhalos/dusthalos/prep_cat_gswlc/overlap_plots'

cat1_basename = 'DoubleMasked_GSWLC-X2_in_SDSS_randoms.fits'
cat2_basename = 'rand_sdss_bg2_JOINED_catalog.fits'
cat1_label = 'GWSLC randoms'
cat2_label = 'SDSS randoms'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(
    cat1_name=cat1, 
    cat2_name=cat2,
    subsample_cat1=True,
    subsample_cat2=True, 
    subsample_size=3e5
)

outname = os.path.join(out_dir, 'rand_gwslc_x_sdss_rectilin.png')
op.make_plot(
    outname=outname, label1=cat1_label, label2=cat2_label,
    coordframe1='icrs', ra_key1='ra', dec_key1='dec', 
    coordframe2='icrs', ra_key2='ra_rand', dec_key2='dec_rand',
)

outname = os.path.join(out_dir, 'rand_gwslc_x_sdss_aitoff.png')
op.make_plot(
    outname=outname, label1=cat1_label, label2=cat2_label,
    coordframe1='icrs', ra_key1='ra', dec_key1='dec', 
    coordframe2='icrs', ra_key2='ra_rand', dec_key2='dec_rand',
    projection='aitoff', central_longitude=180
)

# Save memory!
del(op)
