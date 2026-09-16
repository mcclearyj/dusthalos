from src.plotter import OverlapPlotter
import os

cat_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_csfd'
out_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_csfd/overlap_plots'

if not os.path.exists(out_dir):
    try:
        os.makedirs(out_dir)
    except:
        raise FileNotFoundError(f"Could not create output directory: {out_dir}")

cat1_basename = 'DoubleMasked_wiseScosPhotoz160708_zlt0.2_bCal_gt_16.fits'
cat2_basename = 'redmagic_hidens_z_gt_0.5_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'WISExSCOS galaxies'
cat2_label = 'redMaGiC hi-dens galaxies'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

op = OverlapPlotter(
    cat1_name=cat1, cat2_name=cat2, 
    subsample_cat1=True, subsample_cat2=True, subsample_size=5e5
)

outname = os.path.join(out_dir, 'scos_rmz_hidens_overlap_rectilin.png')
op.make_plot(
    outname=outname, coordframe1='galactic', ra_key1='l', dec_key1='b', 
    coordframe2='icrs', ra_key2='ra_redmagic_hidens_z_gt_0.5', dec_key2='dec_redmagic_hidens_z_gt_0.5',
    label1=cat1_label, label2=cat2_label
)

outname = os.path.join(out_dir, 'scos_rmz_hidens_overlap_aitoff.png')
op.make_plot(
    outname=outname, coordframe1='galactic', ra_key1='l', dec_key1='b', 
    coordframe2='icrs', ra_key2='ra_redmagic_hidens_z_gt_0.5', dec_key2='dec_redmagic_hidens_z_gt_0.5',
    label1=cat1_label, label2=cat2_label, projection='aitoff'
)

# Save memory!
del(op)


###
### Also do it for randoms!
###

cat1_basename = 'DoubleMasked_WISExSCOS_randoms.fits'
cat2_basename = 'redmagic_hidens_randoms_z_gt_0.5_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'WISExSCOS randoms'
cat2_label = 'redMaGiC hi-dens randoms'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

# Randoms are large, so just plot a subsample
op = OverlapPlotter(
    cat1_name=cat1, cat2_name=cat2, 
    subsample_cat1=True, subsample_cat2=True, subsample_size=6e5
)

outname = os.path.join(out_dir, 'rand_scos_rmz_hidens_overlap_rectilin.png')
op.make_plot(
    outname=outname, coordframe1='galactic', ra_key1='l', dec_key1='b',
    coordframe2='icrs', ra_key2='ra', dec_key2='dec', 
    label1=cat1_label, label2=cat2_label
)

outname = os.path.join(out_dir, 'rand_scos_rmz_hidens_overlap_aitoff.png')
op.make_plot(
    outname=outname, coordframe1='galactic', ra_key1='l', dec_key1='b',
    coordframe2='icrs', ra_key2='ra', dec_key2='dec', 
    label1=cat1_label, label2=cat2_label, projection='aitoff'
)

del(op)
