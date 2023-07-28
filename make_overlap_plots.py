from src.plotter import OverlapPlotter

op = OverlapPlotter(cat1_name='../catalogs/WISExSCOS_fg_SCOSmask_cat.fits', cat2_name='../catalogs/wiseSCOS_randoms.fits')

op.make_plot(outname='make_overlap_rectilin.pdf')

op.make_plot(outname='testing_aitoff.png', projection='aitoff')
