from astropy.io import fits
import numpy as np
import os 

fits_table_path = "/n23data1/mccleary/dustyhalos/catalogs/prep_cat_sdss"
fits_table_filename = "rand_sdss_bg2_JOINED_catalog_stacked.fits"
fits_table_file = os.path.join(fits_table_path, fits_table_filename)

with fits.open(fits_table_file) as hdul:
    data = hdul[1].data
    mask = (data['r'] < 22) & (data['r'] > 16)
    newdata = data[mask]
    hdu = fits.BinTableHDU(data=newdata)
    hdu.writeto(
        os.path.join(
            fits_table_path,
            "rand_sdss_bg2_JOINED_catalog_stacked_r_lt_22.fits"
        )
    )
    
# Just to make things easier
datalen = len(data)
hdul.close()
    
# Now make a subsampled catalog of the same length! 
rng = np.random.default_rng()
indices = rng.choice(datalen, size=np.count_nonzero(mask), replace=False)
newdata = data[indices]
hdu = fits.BinTableHDU(data=newdata)
hdu.writeto(
    os.path.join(
        fits_table_path,
        "rand_sdss_bg2_JOINED_catalog_stacked_subsampled.fits"
    )
)

###
### Repeat for galaxies!
###

fits_table_filename = "DoubleMasked_sdss_bg_photoz2.fits"
fits_table_file = os.path.join(fits_table_path, fits_table_filename)

with fits.open(fits_table_file) as hdul:
    data = hdul[1].data
    mask = (data['r'] < 22) & (data['r'] > 16)
    newdata = data[mask]
    hdu = fits.BinTableHDU(data=newdata)
    hdu.writeto(
        os.path.join(
            fits_table_path,
            "DoubleMasked_sdss_bg_photoz2_r_lt_22.fits"
        )
    )

# Just to make things easier                                                                                                  
datalen = len(data)
hdul.close()

# Now make a subsampled catalog of the same length!                                                                           
rng = np.random.default_rng()
indices = rng.choice(datalen, size=np.count_nonzero(mask), replace=False)
newdata = data[indices]
hdu = fits.BinTableHDU(data=newdata)
hdu.writeto(
    os.path.join(
        fits_table_path,
        "DoubleMasked_sdss_bg_photoz2_subsampled.fits"
    )
)
