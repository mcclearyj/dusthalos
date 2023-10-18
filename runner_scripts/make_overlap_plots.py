from src.plotter import OverlapPlotter
import os

#cat_dir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
#outdir = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
cat_dir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_output'                                                                                                      
outdir = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_output'                                                                                                      

cat1_basename = 'Masked_wiseScosPhotoz160708_redshift_cut.fits'
cat2_basename = 'DoubleMasked_Masked_wiseScosPhotoz160708_redshift_cut.fits'
cat3_basename = 'redmagic_hiz_y3_GOLD_JOINED_catalog.fits'

cat1_label = 'WISExSuperCOSMOS'
cat2_label = 'Masked WISExSuperCOSMOS'
cat3_label = 'redMaGiC hi-lum hi-z'

cat1 = os.path.join(cat_dir, cat1_basename)
cat2 = os.path.join(cat_dir, cat2_basename)
cat3 = os.path.join(cat_dir, cat3_basename)

op = OverlapPlotter(cat1_name=cat3, cat2_name=cat2)
op2 = OverlapPlotter(cat1_name=cat1, cat2_name=cat2)

outname = os.path.join(outdir, 'full_scos_hilum_hiz_overlap_rectilin.png')
op.make_plot(outname=outname, 
             ra_tag2='ra_redmagic_hiz', dec_tag2='dec_redmagic_hiz', 
             coordframe2='icrs', label1=cat1_label, label2=cat3_label)

outname = os.path.join(outdir, 'full_scos_hilum_hiz_overlap_mollweide.png')
op.make_plot(outname=outname, projection='mollweide', 
             ra_tag2='ra_redmagic_hiz', dec_tag2='dec_redmagic_hiz', 
             coordframe2='icrs', label1=cat1_label, label2=cat3_label)


outname2 = os.path.join(outdir, 'masked_scos_hilum_hiz_overlap_mollweide.png')
op2.make_plot(outname=outname2, coordframe2='icrs',
              ra_tag2='ra_redmagic_hidens', dec_tag2='dec_redmagic_hidens', 
              label1=cat1_label, label2=cat3_label)

outname2 = os.path.join(outdir, 'masked_scos_hilum_hiz_overlap_mollweide.png')
op2.make_plot(outname=outname2, projection='mollweide',
             ra_tag2='ra_redmagic_hiz', dec_tag2='dec_redmagic_hiz',
             coordframe2='icrs', label1=cat2_label, label2=cat3_label)

