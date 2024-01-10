import numpy as np
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from astropy.table import Table
from scipy import stats
import numpy as np
import pdb

def set_rc_params(fontsize=None):

    print("Setting Matplotlib RC parameters...")

    if fontsize is None:
        fontsize=16
    else:
        fontsize=int(fontsize)

    rc('font',**{'family':'serif'})
    rc('text', usetex=True)

    #plt.rcParams.update({'figure.facecolor':'w'})
    plt.rcParams.update({'axes.linewidth': 1.1})
    plt.rcParams.update({'xtick.labelsize': fontsize})
    plt.rcParams.update({'ytick.labelsize': fontsize})
    plt.rcParams.update({'xtick.major.size': 8})
    plt.rcParams.update({'xtick.major.width': 1.1})
    plt.rcParams.update({'xtick.minor.visible': True})
    plt.rcParams.update({'xtick.minor.width': 1.})
    plt.rcParams.update({'xtick.minor.size': 6})
    plt.rcParams.update({'xtick.direction': 'out'})
    plt.rcParams.update({'ytick.major.width': 1.1})
    plt.rcParams.update({'ytick.major.size': 8})
    plt.rcParams.update({'ytick.minor.visible': True})
    plt.rcParams.update({'ytick.minor.width': 1.})
    plt.rcParams.update({'ytick.minor.size':6})
    plt.rcParams.update({'ytick.direction':'out'})
    plt.rcParams.update({'axes.labelsize': fontsize})
    plt.rcParams.update({'axes.titlesize': fontsize})
    plt.rcParams.update({'legend.fontsize': int(fontsize-2)})

    return

def remove_outliers(data, bands):
    # All bands are good <3
    catlen = len(data)
    wg = np.full(catlen, True)

    # Pick out only the good entries!
    for band in bands:
        band_mag = np.ma.getdata(data[band])
        band_bool = (band_mag > -9999) & (band_mag != np.nan) #& (band_mag < 27)
        wg *= band_bool

    # percent of galaxies that failed
    pfail = 100-(np.count_nonzero(wg) / catlen * 100)

    # How many failed?
    print(f'RedshiftCalc: {np.count_nonzero(wg)}/{catlen} galaxies',
          f'({100-pfail:.1f}%) have good photometry')
    print(f'RedshiftCalc: removing {pfail:.1f}% of galaxies from data')
    print('')

    return wg

def bin_the_redshifts(rmz_rand, n_zbins=10):
    ## Create colors
    gmr = rmz_rand['mof_cm_mag_corrected_g'] - rmz_rand['mof_cm_mag_corrected_r']
    imz = rmz_rand['mof_cm_mag_corrected_i'] - rmz_rand['mof_cm_mag_corrected_z']

    ## This is kinda like the algorithm I use elsewhere
    zhist, zbins = np.histogram(rmz_rand['z'], bins=n_zbins)
    
    # Hold redshift bins
    zbin_values = np.digitize(rmz_rand['z'], bins=zbins, right=False)
    
    # Set the number of bins
    bin_numbers = np.unique(zbin_values)
    
    # Hold mean color values & redshift
    gmr_values = []; imz_values = []; z_bin_median = []

    # Loop over all redshift bins specified in bin_numbers
    for zb in bin_numbers:
        
        # select only galaxies in the current redshift bin
        c_slice = (zbin_values == zb)
        this_gmr_bin = gmr[c_slice]
        this_imz_bin = imz[c_slice]
        this_z_bin = rmz_rand['z'][c_slice]

        gmr_values.append((np.mean(this_gmr_bin),
                           np.std(this_gmr_bin)
                           ))
        imz_values.append((np.mean(this_imz_bin),
                           np.std(this_imz_bin)
                           ))
        z_bin_median.append(np.mean(this_z_bin))

    tab = Table([z_bin_median,
                 np.array(gmr_values)[:,0], np.array(gmr_values)[:,1],
                 np.array(imz_values)[:,0], np.array(imz_values)[:,1]],
                 names=['z_median', 'gmr_median', 'gmr_std',
                 'imz_median', 'imz_std'])

    tab.write('binned_table.txt', format='ascii.csv', overwrite=True)

    return tab

def main():
    #catalog_path = '/work/mccleary_group/dusty_halos/catalogs/prep_cat_output'
    catalog_path = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
    rmz_cat = Table.read(os.path.join(catalog_path,
                     'z_lt_045_redmagic_hidens_y3_GOLD_JOINED_catalog.fits'),
                     memmap=True)
    rmz_rand_f = fits.open(os.path.join(catalog_path,
                         'z_lt_045_redmagic_hidens_randoms_y3_GOLD_JOINED_catalog.fits'))
    
    rng = np.random.default_rng()
    rint = rng.integers(0, high=len(rmz_rand_f[1].data), size=int(2e6), dtype=np.int64)
    rmz_rand=rmz_rand_f[1].data[rint]

    bands = ['mof_cm_mag_corrected_g', 'mof_cm_mag_corrected_r',
             'mof_cm_mag_corrected_i', 'mof_cm_mag_corrected_z']

    wg_rmz_cat = remove_outliers(rmz_cat, bands)
    rmz_cat = rmz_cat[wg_rmz_cat]
    wg_rmz_rand = remove_outliers(rmz_rand, bands)
    rmz_rand = rmz_rand[wg_rmz_rand]

    binned_randoms = bin_the_redshifts(rmz_rand)

    fig, axs = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)

    axs[0].plot(rmz_cat['zredmagic'],
                (rmz_cat['mof_cm_mag_corrected_g'] - rmz_cat['mof_cm_mag_corrected_r']),
                '.', markersize=0.5, color='xkcd:deep red', label='hidens $z < 0.45$')

    axs[0].errorbar(binned_randoms['z_median'], binned_randoms['gmr_median'],
                    yerr=binned_randoms['gmr_std'], fmt='o-',
                    color='xkcd:tomato red', label='hidens $z < 0.45$ randoms',
                    capsize=5)
    
    axs[0].set_ylim(-1,4)
    axs[0].legend(loc='upper left')
    axs[0].set_xlabel('Redshift')
    axs[0].set_ylabel(r'$g$ - $r$')

    axs[1].plot(rmz_cat['zredmagic'],
                (rmz_cat['mof_cm_mag_corrected_i'] - rmz_cat['mof_cm_mag_corrected_z']),
                '.', markersize=0.5, color='xkcd:deep red', label='hidens $z < 0.45$')

    axs[1].errorbar(binned_randoms['z_median'], binned_randoms['imz_median'],
                    yerr=binned_randoms['imz_std'], fmt='o-',
                    color='xkcd:tomato red', label='hidens $z < 0.45$ randoms',
                    capsize=5)

    axs[1].set_ylim(-1,2)
    axs[1].legend(loc='upper left')
    axs[1].set_xlabel('Redshift')
    axs[1].set_ylabel(r'$i$ - $z$')

    fig.suptitle('Redmagic z $<$ 0.45 resampled') 
    fig.savefig('color_redshift_redmagic_z_lt_045_hidens_resamp.png')

    return 0

if __name__ == '__main__':

    rc = main()

    if rc == 0:
        print("color redshift plots have been successfully completed")
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
