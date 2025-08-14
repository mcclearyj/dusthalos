from src.plotter import OverlapPlotter
import os

cat_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_csfd'
out_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_csfd/overlap_plots'
if not os.path.exists(out_dir):
    try:
        os.makedirs(out_dir)
    except:
        raise FileNotFoundError(f"Could not create output directory: {out_dir}")

cat1_basename = 'DoubleMasked_gaia_stars_southern_subsample1.fits'
cat2_basename = 'redmagic_hiz_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'Gaia stars'
cat2_label = 'redMaGiC hi-z galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(
    cat1_name=cat1, cat2_name=cat2, 
    subsample_cat1=False, subsample_cat2=True, subsample_size=6e5
)

outname = os.path.join(out_dir, 'gaia_rmz_hiz_overlap_rectilin.png')
op.make_plot(
    outname=outname, coordframe1='icrs', ra_key1='ra', dec_key1='dec', 
    coordframe2='icrs', ra_key2='ra_redmagic_hiz', dec_key2='dec_redmagic_hiz',
    label1=cat1_label, label2=cat2_label
)

outname = os.path.join(out_dir, 'gaia_rmz_hiz_overlap_aitoff.png')
op.make_plot(
    outname=outname, coordframe1='icrs', ra_key1='ra', dec_key1='dec', 
    coordframe2='icrs', ra_key2='ra_redmagic_hiz', dec_key2='dec_redmagic_hiz',
    label1=cat1_label, label2=cat2_label, projection='aitoff'
)

# Save memory!
del(op)


###
### Also do it for randoms!
###

cat1_basename = 'DoubleMasked_gaia_stars_southern_subsample2.fits'
cat2_basename = 'redmagic_hiz_randoms_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'Gaia randoms'
cat2_label = 'redMaGiC hi-z randoms'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

# Randoms are large, so just plot a subsample
op = OverlapPlotter(
    cat1_name=cat1, cat2_name=cat2, 
    subsample_cat1=False, subsample_cat2=True, subsample_size=6e5
)

outname = os.path.join(out_dir, 'rand_gaia_rmz_hiz_overlap_rectilin.png')
op.make_plot(
    outname=outname, coordframe1='icrs', ra_key1='ra', dec_key1='dec',
    coordframe2='icrs', ra_key2='ra', dec_key2='dec', 
    label1=cat1_label, label2=cat2_label
)

outname = os.path.join(out_dir, 'rand_gaia_rmz_hiz_overlap_aitoff.png')
op.make_plot(
    outname=outname, coordframe1='icrs', ra_key1='ra', dec_key1='dec',
    coordframe2='icrs', ra_key2='ra', dec_key2='dec', 
    label1=cat1_label, label2=cat2_label, projection='aitoff'
)

del(op)
