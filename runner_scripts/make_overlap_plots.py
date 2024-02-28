from src.plotter import OverlapPlotter
import os

cat_dir = '/work/mccleary_group/dusty_halos/catalogs/'
out_dir = '/work/mccleary_group/dusty_halos/catalogs/'

cat1_basename = 'GSWLC-X2_in_SDSS.fits'
cat2_basename = 'prep_cat_output/redmagic_hiz_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'GSWLC-X2 in SDSS'
cat2_label = 'redMaGiC hi-z'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)

# Make another one
op = OverlapPlotter(cat1_name=cat1, cat2_name=cat2)

outname = os.path.join(out_dir, 'GSWLC-X2_redmagic_overlap_rectilin.pdf')
op.make_plot(
    outname=outname, label1=cat1_label, coordframe1='icrs',
    ra_tag1='RA', dec_tag1='Dec', ra_tag2='ra_redmagic_hiz', 
    dec_tag2='dec_redmagic_hiz', label2=cat2_label, coordframe2='icrs'
)

outname = os.path.join(out_dir, 'GSWLC-X2_redmagic_overlap_aitoff.pdf')
op.make_plot(
    outname=outname, label1=cat1_label, projection='aitoff', 
    coordframe1='icrs', ra_tag1='RA', dec_tag1='Dec', ra_tag2='ra_redmagic_hiz', 
    dec_tag2='dec_redmagic_hiz', label2=cat2_label, coordframe2='icrs'
)
             
             
