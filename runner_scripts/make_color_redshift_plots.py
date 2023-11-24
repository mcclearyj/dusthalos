import numpy as np
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from astropy.table import Table
from scipy import stats


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
        band_bool = (band_mag > -9999) & (band_mag != np.nan) & (band_mag < 27)
        wg *= band_bool

    # percent of galaxies that failed
    pfail = 100-(np.count_nonzero(wg) / catlen * 100)

    # How many failed?
    print(f'RedshiftCalc: {np.count_nonzero(wg)}/{catlen} galaxies',
          f'({100-pfail:.1f}%) have good photometry')
    print(f'RedshiftCalc: removing {pfail:.1f}% of galaxies from data')
    print('')

    return wg


def main():
    catalog_path = '/Users/j.mccleary/Research/dusty_halos/catalogs/prep_cat_output'
    hiz = Table.read(os.path.join(catalog_path, 'redmagic_hiz_y3_GOLD_JOINED_catalog.fits'), memmap=True)
    #hidens = Table.read(os.path.join(catalog_path,'redmagic_hidens_y3_GOLD_JOINED_catalog.fits'))
    hidens = Table.read(os.path.join(catalog_path,'redmagic_hiz_randoms_y3_GOLD_JOINED_catalog.fits'), memmap=True)

    bands = ['mof_cm_mag_corrected_g', 'mof_cm_mag_corrected_r',
        'mof_cm_mag_corrected_i', 'mof_cm_mag_corrected_z']


    wg_hiz = remove_outliers(hiz, bands)
    hiz = hiz[wg_hiz]
    wg_hidens = remove_outliers(hidens, bands)
    hidens = hidens[wg_hidens]


    fig, axs = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)
    axs[0].plot(hidens['zredmagic'], 
            (hidens['mof_cm_mag_corrected_g'] - hidens['mof_cm_mag_corrected_r']), 
            '.', markersize=0.5, color='xkcd:tomato red', label='hilum_hiz rand')
    axs[0].plot(hiz['zredmagic'], 
            (hiz['mof_cm_mag_corrected_g'] - hiz['mof_cm_mag_corrected_r']), 
             '.', markersize=0.5, color='xkcd:deep red', label='hilum_hiz')
    axs[0].set_ylim(-1,6)
    axs[0].legend(markerscale=15, loc='upper left')
    axs[0].set_xlabel('Redshift')
    axs[0].set_ylabel(r'$g$ - $r$')

    axs[1].plot(hidens['zredmagic'], 
                (hidens['mof_cm_mag_corrected_i'] - hidens['mof_cm_mag_corrected_z']), 
                '.', markersize=0.5, color='xkcd:tomato red', label='hilum_hiz rand')
    axs[1].plot(hiz['zredmagic'], 
                (hiz['mof_cm_mag_corrected_i'] - hiz['mof_cm_mag_corrected_z']), 
                 '.', markersize=0.5, color='xkcd:deep red', label='hilum_hiz')
    axs[1].set_ylim(-0.5,1.5)
    axs[1].legend(markerscale=15, loc='upper left')
    axs[1].set_xlabel('Redshift')
    axs[1].set_ylabel(r'$i$ - $z$')

    fig.savefig('color_redshift_redmagic_hiz_rand.png')

    return 0 

if __name__ == '__main__':
    
    rc = main()

    if rc == 0:
        print("color redshift plots have been successfully completed")
    else:
        print(f'get_jwst_psfs.py has failed w/ rc={rc}')
    
