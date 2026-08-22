import healpy as hp
import matplotlib.pyplot as plt
from astropy.table import Table
import os
from src.hpmask import HpMask
import numpy as np

output_dir = "/projects/mccleary_group/hsia.i/dusthalos_output/"
output_name = "demeaned_sdss.fits"
input_dir = "/projects/mccleary_group/dusty_halos/dusthalos_emh/output"
catname = "sdss_csfd/dustcorrel_sdss_bg_photoz2_treecorrcat.fits"

sdss = Table.read(os.path.join(input_dir,catname), memmap=True, format="fits")

ra_sdss = sdss["ra"].data
dec_sdss = sdss["dec"].data
av_sdss = sdss["k"].data

cat_pix = HpMask.coords_to_healpixels(lon=ra_sdss, lat=dec_sdss, frame="icrs", nside=256)
map_file = HpMask(os.path.join(input_dir, "sdss_mean_av_map.fits"), coordframe="icrs")

demeaned_av_sdss = av_sdss - map_file.mask[cat_pix]

out = Table()
out["ra"] = ra_sdss
out["dec"] = dec_sdss
out["av"] = av_sdss
out["demeaned_av"] = demeaned_av_sdss
out.write(os.path.join(output_dir, output_name), overwrite=True)
print("Demeaning complete!")
