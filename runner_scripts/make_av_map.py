from src.plotter import AvMapPlotter
import os

cat_dir = '/projects/mccleary_group/dusty_halos/dusthalos_emh/output'
out_dir = '/projects/mccleary_group/dusty_halos/dusthalos_emh/output'

## Start with full SDSS treecorrcat
cat1_name = os.path.join(
    cat_dir, 
    'sdss_csfd/dustcorrel_sdss_bg_photoz2_treecorrcat.fits'
)
outname1 = os.path.join(cat_dir, 'sdss_csfd/sdss_bg_photoz2_av_map.pdf')
plot_title1 = r"Mean TreeCorr $A_V$ of SDSS Background Galaxies"

## Then try the WISExSCOS x redMaGiC
cat2_name = os.path.join(
    cat_dir, 
    'redmagic_hidens_csfd/dustcorrel_demeancov_fixed_est_treecorrcat.fits'
)
outname2 = os.path.join(cat_dir, 'redmagic_hidens_csfd/hidens_csfd_av_map.pdf')
plot_title2 = r"Mean TreeCorr $A_V$ of WISExSCOS x redMaGiC Hi-Dens Background Galaxies"

# Make map 1
av_map = AvMapPlotter(cat_name=cat1_name)
av_map.make_Av_map(out_name=outname1, plot_title=plot_title1)

# Make map 2
av_map = AvMapPlotter(cat_name=cat2_name, nside=512)
av_map.make_Av_map(out_name=outname2, plot_title=plot_title2)



