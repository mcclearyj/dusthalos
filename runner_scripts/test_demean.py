import healpy as hp
import matplotlib.pyplot as plt
from astropy.table import Table
import os
from src.hpmask import HpMask
import numpy as np

input_dir = "/projects/mccleary_group/dusty_halos/dusthalos_emh/output"

sdss = Table.read(
    os.path.join(input_dir, "sdss_csfd/dustcorrel_sdss_bg_photoz2_treecorrcat.fits"),
    memmap=True, format="fits")

ra_sdss = sdss["ra"].data
dec_sdss = sdss["dec"].data
av_sdss = sdss["k"].data

cat_pix = HpMask.coords_to_healpixels(lon=ra_sdss, lat=dec_sdss, frame="icrs", nside=256)
map_file = HpMask(os.path.join(input_dir, "sdss_csfd/sdss_mean_av_map.fits"), coordframe="icrs")

demeaned_av_sdss = av_sdss - map_file.mask[cat_pix]

out = Table()
out["ra"] = ra_sdss
out["dec"] = dec_sdss
out["av"] = av_sdss
out["demeaned_av"] = demeaned_av_sdss
out.write("/projects/mccleary_group/hsia.i/dusthalos_output/demeaned_sdss.fits", overwrite=True)
print("Demeaning complete!")
