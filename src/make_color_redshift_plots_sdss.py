import numpy as np
from matplotlib import rc,rcParams
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from astropy.table import Table
from scipy import stats
import numpy as np
import pdb
import seaborn as sns

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
    plt.rcParams.update({'legend.fontsize': int(fontsize)})

    return

def remove_outliers(data, bands):
    # All bands are good <3
    catlen = len(data)
    wg = np.full(catlen, True)

    # Pick out only the good entries!
    for band in bands:
        band_mag = np.ma.getdata(data[band])
        band_bool = (np.abs(band_mag) < 99) & (band_mag != np.nan) & (band_mag < 27)
        wg *= band_bool

    # percent of galaxies that failed
    pfail = 100-(np.count_nonzero(wg) / catlen * 100)

    # How many failed?
    print(f'RedshiftCalc: {np.count_nonzero(wg)}/{catlen} galaxies',
          f'({100-pfail:.1f}%) have good photometry')
    print(f'RedshiftCalc: removing {pfail:.1f}% of galaxies from data')
    print('')

    return wg

def bin_the_redshifts(rand_cat, bands, z2_colname, n_zbins=17):
    ## Create colors
    umg = rand_cat[bands[0]] - rand_cat[bands[1]]
    gmr = rand_cat[bands[1]] - rand_cat[bands[2]]
    imz = rand_cat[bands[3]] - rand_cat[bands[4]]

    ## This is kinda like the algorithm I use elsewhere
    zhist, zbins = np.histogram(rand_cat[z2_colname], bins=n_zbins)

    # Hold redshift bins
    zbin_values = np.digitize(rand_cat[z2_colname], bins=zbins, right=False)

    # Set the number of bins
    bin_numbers = np.unique(zbin_values)

    # Hold mean color values & redshift
    umg_values = []; gmr_values = []; imz_values = []; z_bin_median = []

    # Loop over all redshift bins specified in bin_numbers
    for zb in bin_numbers:

        # select only galaxies in the current redshift bin
        c_slice = (zbin_values == zb)
        sq_ngals = np.sqrt(np.count_nonzero(c_slice))
        this_umg_bin = umg[c_slice]
        this_gmr_bin = gmr[c_slice]
        this_imz_bin = imz[c_slice]
        this_z_bin = rand_cat[z2_colname][c_slice]

        umg_values.append((np.median(this_umg_bin), np.std(this_umg_bin)))
        gmr_values.append((np.median(this_gmr_bin), np.std(this_gmr_bin)))
        imz_values.append((np.median(this_imz_bin), np.std(this_imz_bin)))
        z_bin_median.append(np.median(this_z_bin))

    tab = Table(
        [
            z_bin_median, 
            np.array(umg_values)[:,0], np.array(umg_values)[:,1],
            np.array(gmr_values)[:,0], np.array(gmr_values)[:,1],
            np.array(imz_values)[:,0], np.array(imz_values)[:,1]
        ],
        names=[
            'z_median', 'umg_median', 'umg_std', 
            'gmr_median', 'gmr_std', 'imz_median', 'imz_std'
        ]
    )

    tab.write('binned_color_table.txt', format='ascii.csv', overwrite=True)

    return tab

def main():
    # Make plots look pretty
    set_rc_params(fontsize=16)

    # Catalog path
    catalog_path = '/n23data1/mccleary/dustyhalos/catalogs/prep_cat_sdss'


    # Read in galaxies
    gals_cat = Table.read(
        os.path.join(catalog_path,
        'DoubleMasked_sdss_bg_photoz2_r_lt_22.fits'), memmap=True
    )

    # Read in random catalog
    rand_cat_f = fits.open(os.path.join(catalog_path,
        'rand_sdss_bg2_JOINED_catalog_stacked_r_lt_22.fits'), memmap=True
    )


    # Status update; should be logger but w/e
    print("Catalogs read in...")
    
    # In case catalogs are super-super long, pick subset for plotting
    if len(gals_cat) > int(6e5):
        # Status update; should be logger but w/e
        print("Subsampling galaxy catalog...")
        rng = np.random.default_rng()
        rint = rng.integers(
            0, high=len(gals_cat),
            size=int(6e5), dtype=np.int64
        )
        gals_cat = gals_cat[rint]
    
    if len(rand_cat_f[1].data) > int(6e5):
        print("Subsampling random catalog...")
        rng = np.random.default_rng()
        rint = rng.integers(0, high=len(rand_cat_f[1].data),
                size=int(6e5), dtype=np.int64)
        rand_cat=rand_cat_f[1].data[rint]

    else:
        rand_cat=rand_cat_f[1].data

    bands = [
        'u_corr_csfd', 'g_corr_csfd', 'r_corr_csfd', 
        'i_corr_csfd', 'z_corr_csfd'
    ]

    z1_colname = 'redshift'
    z2_colname = 'redshift_rand'

    wg_gals_cat = remove_outliers(gals_cat, bands)
    gals_cat = gals_cat[wg_gals_cat]
    wg_rand_cat = remove_outliers(rand_cat, bands)
    rand_cat = rand_cat[wg_rand_cat]

    print("Binning random catalog...")
    binned_randoms = bin_the_redshifts(rand_cat, bands=bands, z2_colname=z2_colname)
    
    print("Making color-redshift scatterplot...")
    fig, axs = plt.subplots(3, 1, figsize=(6.5, 18), tight_layout=True)

    axs[0].plot(gals_cat[z1_colname],
                (gals_cat[bands[0]] - gals_cat[bands[1]]),
                '.', markersize=0.04, color='xkcd:deep red', label='Galaxies')

    axs[0].errorbar(
        binned_randoms['z_median'], binned_randoms['umg_median'],
        yerr=binned_randoms['umg_std'], fmt='o-', color='xkcd:tomato red',
        label='Randoms', capsize=5
    )

    axs[0].set_ylim(-3, 5)
    lgnd = axs[0].legend(loc='upper right', scatterpoints=1)
    # Makes galaxy points the same size as errplot points! 
    lgnd.legend_handles[0].set_markersize(10)  
    axs[0].set_xlabel('Redshift')
    axs[0].set_ylabel(r'$u - g$')


    axs[1].plot(gals_cat[z1_colname],
                (gals_cat[bands[1]] - gals_cat[bands[2]]),
                '.', markersize=0.04, color='xkcd:deep red', label='Galaxies')

    axs[1].errorbar(
        binned_randoms['z_median'], binned_randoms['gmr_median'],
        yerr=binned_randoms['gmr_std'], fmt='o-', color='xkcd:tomato red',
        label='Randoms', capsize=5
    )

    axs[1].set_ylim(-0.6, 3.5)
    lgnd = axs[1].legend(loc='upper right', scatterpoints=1)
    # Makes galaxy points the same size as errplot points! 
    lgnd.legend_handles[0].set_markersize(10)  
    axs[1].set_xlabel('Redshift')
    axs[1].set_ylabel(r'$g - r$')

    axs[2].plot(
        gals_cat[z1_colname], (gals_cat[bands[3]] - gals_cat[bands[4]]),
        '.', markersize=0.04, color='xkcd:deep red', label='Galaxies'
    )

    axs[2].errorbar(
        binned_randoms['z_median'], binned_randoms['imz_median'],
        yerr=binned_randoms['imz_std'], fmt='o-', color='xkcd:tomato red',
        label='Randoms', capsize=5
    )

    axs[2].set_ylim(-2.5, 2.5)
    lgnd = axs[2].legend(loc='upper right')
    lgnd.legend_handles[0].set_markersize(10)
    axs[2].set_xlabel('Redshift')
    axs[2].set_ylabel(r'$i - z$')

    fig.suptitle(r'SDSS Photoz2 $r < 22$ sample', fontsize=14)
    fig.savefig('color_redshift_sdss_photoz2_r_lt_22_resamp.png')
    print("Scatterplot made!")
    
    ## Try a Seaborn plot
    
    gals_pd = gals_cat.to_pandas()
    sns.set_style("whitegrid")
    print("Making contour plot...")
    fig, axs = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)
    sns.kdeplot(
        x=gals_pd[z1_colname],
        y=(gals_pd[bands[0]] - gals_pd[bands[1]]),
        fill=True, 
        cmap="Blues",
        levels=10,  
        thresh=0.05, 
        clip=[[0.5, 1], [-2.5, 4.5]],
        ax=axs[0]
    )
    axs[0].set_xlabel(z1_colname)
    axs[0].set_ylabel(f"{bands[0]} - {bands[1]}")
    axs[0].set_ylim(-3, 5)
    
    sns.kdeplot(
        x=gals_pd[z1_colname],
        y=(gals_pd[bands[2]] - gals_pd[bands[3]]),
        fill=True, 
        cmap="Blues",
        levels=10, 
        thresh=0.05, 
        ax=axs[1], 
    )
    axs[1].set_xlabel(z1_colname)
    axs[1].set_ylabel(f"{bands[2]} - {bands[3]}")
    axs[1].set_ylim(-1.5, 2.5)

    fig.suptitle(r'SDSS Photoz2 $r < 22$ sample')
    fig.savefig('color_redshift_sdss_photoz2_r_lt_22_contour.png')
    print("Contour plot made!")
    
    return 0

if __name__ == '__main__':

    rc = main()

    if rc == 0:
        print("color redshift plotting has successfully completed")
    else:
        print(f'color redshift plotting has failed w/ rc={rc}')
