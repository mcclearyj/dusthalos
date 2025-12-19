import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
import os

# === CONFIGURATION ===
nside = 64
npix = hp.nside2npix(nside)

# === BUILD HEALPIX MAPS ===
def coords_to_healpixels(lon, lat, frame):
    # Build SkyCoord object from input latitude and longitude
    if frame == 'icrs':
        coords = SkyCoord(ra=lon*u.deg, dec=lat*u.deg, frame='icrs')
    elif frame == 'galactic':
        coords_gal = SkyCoord(l=lon*u.deg, b=lat*u.deg, frame='galactic')
        coords = coords_gal.icrs
    # Assign galaxies to HEALPixels for given NSIDE
    theta = np.radians(90.0 - coords.dec.deg)
    phi = np.radians(coords.ra.deg)
    pixels = hp.ang2pix(nside, theta, phi)
    return pixels

def build_mean_map(lon, lat, values, frame='icrs'):
    pix = coords_to_healpixels(lon, lat, frame)
    sums = np.bincount(pix, weights=values, minlength=npix)
    counts = np.bincount(pix, minlength=npix)
    result = np.full(npix, hp.UNSEEN)
    result[counts > 0] = sums[counts > 0] / counts[counts > 0]
    return result

def build_density_map(lon, lat, frame='icrs'):
    pix = coords_to_healpixels(lon, lat, values)
    # Tally up counts
    counts = np.bincount(pix, minlength=npix)
    result = np.full(npix, hp.UNSEEN)
    result[counts > 0] = counts[counts > 0]
    return result

# === LOAD DATA ===
cat_dir = "/Users/j.mccleary/Research/dusty_halos/catalogs"
outputs_dir = "/Users/j.mccleary/Research/dusty_halos/dusthalos_emh/output"

wise = Table.read(
    os.path.join(cat_dir,
    "wiseScosPhotoz160708_zlt0.2_bCal_gt_16.fits"
    )
)
wise_sdss = Table.read(
    os.path.join(cat_dir,
    "prep_cat_sdss/DoubleMasked_wiseScosPhotoz160708_zlt0.15_rCal_gt_17.fits"
    )
)
wise_rm = Table.read(
    os.path.join(cat_dir,
    "prep_cat_csfd/DoubleMasked_wiseScosPhotoz160708_zlt0.2_bCal_gt_16.fits"
    )
)
redmagic = Table.read(
    os.path.join(outputs_dir,
    "redmagic_hidens_csfd/dustcorrel_demeancov_fixed_est_treecorrcat.fits"
    )
)
sdss = Table.read(
    os.path.join(outputs_dir,
    "sdss_csfd/dustcorrel_sdss_bg_photoz2_treecorrcat.fits"
    ), memmap=True, format="fits"
)

sdss_fg = Table.read(
    os.path.join(cat_dir,
    "prep_cat_sdss/DoubleMasked_sdss_fg_photoz.fit"
    ), memmap=True, format="fits"
)

vector = 'k'
ra_sdss = sdss['ra'].data
dec_sdss = sdss['dec'].data
av_sdss = sdss[vector].data

ra_sdss_fg = sdss_fg['ra'].data
dec_sdss_fg = sdss_fg['dec'].data

ra_redmagic = redmagic['ra'].data
dec_redmagic = redmagic['dec'].data
av_redmagic = redmagic['k'].data

l_wise = wise['l'].data
b_wise = wise['b'].data

l_wise_sdss = wise_sdss['l'].data
b_wise_sdss = wise_sdss['b'].data

l_wise_rm = wise_rm['l'].data
b_wise_rm = wise_rm['b'].data

map_sdss_av = build_mean_map(ra_sdss, dec_sdss, av_sdss)
map_wise_sdss_density = build_density_map(l_wise_sdss, b_wise_sdss, frame='galactic')

map_fg_sdss_density = build_density_map(ra_sdss_fg, dec_sdss_fg, frame='icrs')

map_redmagic_av = build_mean_map(ra_redmagic, dec_redmagic, av_redmagic)
map_wise_rm_density = build_density_map(l_wise_rm, b_wise_rm, frame='galactic')

# === MASKING ===
mask_sdss = (map_sdss_av != hp.UNSEEN) & (map_fg_sdss_density != hp.UNSEEN)
map_sdss_av_masked = np.copy(map_sdss_av)
map_sdss_av_masked[~mask_sdss] = 0
map_fg_sdss_density_masked = np.copy(map_fg_sdss_density)
map_fg_sdss_density_masked[~mask_sdss] = 0

mask_redmagic = (map_redmagic_av != hp.UNSEEN) & (map_wise_rm_density != hp.UNSEEN)
map_redmagic_av_masked = np.copy(map_redmagic_av)
map_redmagic_av_masked[~mask_redmagic] = 0
map_wise_rm_density_masked = np.copy(map_wise_rm_density)
map_wise_rm_density_masked[~mask_redmagic] = 0

# === SPHERICAL HARMONICS ===
alm_sdss = hp.map2alm(map_sdss_av_masked, lmax=200)
alm_redmagic = hp.map2alm(map_redmagic_av_masked, lmax=200)
alm_wise_sdss = hp.map2alm(map_fg_sdss_density_masked, lmax=200)
alm_wise_redmagic = hp.map2alm(map_wise_rm_density_masked, lmax=200)

cl_cross_sdss = hp.alm2cl(alm_sdss, alm_wise_sdss)
cl_cross_redmagic = hp.alm2cl(alm_redmagic, alm_wise_redmagic)

# === PLOT ===
ell = np.arange(len(cl_cross_sdss))
plt.figure(figsize=(10, 6))
plt.plot(ell, cl_cross_sdss, label=r'SDSS $A_V \times$  WISE density', color='blue')
plt.plot(ell, cl_cross_redmagic, label=r'redMaGiC $A_V \times$ WISE density', color='green')
plt.xlabel(r'Multipole $\ell$')
plt.ylabel(r'$C_\ell$')
plt.yscale('log')
plt.title(r'Cross-Power Spectra: $A_V \times$ WISE Galaxy Density')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("cross_power_spectra_wise_sdss_rm.pdf")

# === AV AND DENSITY MAP PLOTS ===
fig = plt.figure(figsize=(9, 5.5))
central_longitude = 0
plot_title = r'SDSS A_V'
hp.mollview(
    map_sdss_av, title=plot_title,
    unit='mag', min=np.nanpercentile(map_sdss_av[mask_sdss], 1),
    max=np.nanpercentile(map_sdss_av[mask_sdss], 99),
    cmap='viridis', fig=1, xsize=1100,rot=[central_longitude, 0, 0],
)

fig = plt.figure(figsize=(9, 5.5))
central_longitude = 0
plot_title = r'WISExSCOS x SDSS galaxy density'
hp.mollview(
    map_wise_sdss_density, title=plot_title,
    unit='Number',
    min=np.nanpercentile(map_fg_sdss_density[map_wise_sdss_density!=hp.UNSEEN], 1),
    max=np.nanpercentile(map_fg_sdss_density[map_wise_sdss_density!=hp.UNSEEN], 98),
    cmap='viridis', xsize=1100, rot=[central_longitude, 0, 0],
)

fig = plt.figure(figsize=(9, 5.5))
central_longitude = 0
plot_title = r'SDSS fg galaxy density'
hp.mollview(
    map_fg_sdss_density_masked, title=plot_title,
    unit='Number',
    min=np.nanpercentile(map_fg_sdss_density[mask_sdss], 1),
    max=np.nanpercentile(map_fg_sdss_density[mask_sdss], 98),
    cmap='viridis', fig=5, xsize=1100, rot=[central_longitude, 0, 0],
)

fig = plt.figure(figsize=(9, 5.5))
central_longitude = 180
plot_title = r'WISExSCOS galaxy density'
hp.mollview(
    map_wise_sdss_density, title=plot_title,
    unit='mag',
    min=np.nanpercentile(map_wise_sdss_density[map_wise_sdss_density != hp.UNSEEN], 1),
    max=np.nanpercentile(map_wise_sdss_density[map_wise_sdss_density != hp.UNSEEN], 98),
    cmap='viridis', xsize=1100, rot=[central_longitude, 0, 0]
)
