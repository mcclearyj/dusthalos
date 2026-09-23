import os, sys
import time
import argparse
import treecorr
import pdb
import numpy as np

# Local imports
import src.utils as utils
from src.correlator import Correlator
from src.plotter import DustPlotter


'''
Under construction! I haven't decided whether this should be a class or just
keep it as is.

General algorithm:
    - Create instances of Correlator w/ correlator configs, load cats
    - Calculate anomalous Av/reddening dust with optimal estimator
    - Do correlation
'''

def set_treecorr_threads(correl_config):
    '''
    TreeCorr's pair counting is OpenMP-parallel, but with num_threads unset it
    sizes the thread pool from the machine's core count. Under SLURM that's the
    whole node, so a job holding a fraction of a node ends up oversubscribing
    its cores. Take the allocation instead, unless the config says otherwise.
    '''
    params = correl_config['treecorr_params']

    if params.get('num_threads') is None:
        params['num_threads'] = utils.get_n_cpus()

    print(f"TreeCorr will use num_threads = {params['num_threads']}\n")


def make_names(correl_config):
    outdir = correl_config['output_path']
    base = correl_config['output_basename']

    names = {'dk_outfile': os.path.join(outdir, base + '_raw_signal.txt'),
             'dr_outfile': os.path.join(outdir, base + '_bg_randoms.txt'),
             'fr_outfile': os.path.join(outdir, base + '_fg_randoms.txt'),
             'rr_outfile': os.path.join(outdir, base + '_fgxbg_randoms.txt'),
             'ck_outfile': os.path.join(outdir, base + '_compensated_signal.txt'),
             'cov_output': os.path.join(outdir, base + '_covariance.txt'),
             'fig_output': os.path.join(outdir, base + '_figure.png')
             }

    return utils.AttrDict(names)

def get_mean_fg_redshift(fg, correl_config):
    '''
    Median redshift of the foreground sample, used to convert angular
    separations to physical ones. Saved treecorr catalogs don't carry a
    redshift column, so fall back to arcminute scaling in that case.
    '''
    if fg.from_file == True:
        print("Foreground loaded from saved treecorr catalog, no redshifts")
        print("Forcing arcminute plot scaling")
        correl_config['use_kpc'] = False
        return 0

    try:
        return np.median(fg.Catalog.data[fg.cat_config['z_key']])
    except KeyError:
        print("No redshift column found, setting mean fg redshift to 0")
        print("Forcing arcminute plot scaling")
        correl_config['use_kpc'] = False
        return 0


def check_patch_consistency(correl_config, **catalogs):
    '''
    Jackknife and sample covariances require every catalog to share the same
    patches. That's automatic when the patches are assigned from a common set
    of patch centers, but catalogs read from file bring their own patch column,
    so two files written by different runs would silently disagree.
    '''
    var_method = correl_config['treecorr_params'].get('var_method')

    if var_method not in ['jackknife', 'sample']:
        return

    npatches = {name: cat.treecorrCatalog.npatch
                    for name, cat in catalogs.items()}

    if len(set(npatches.values())) > 1:
        raise ValueError(
            f"var_method '{var_method}' requires consistent patches across " + \
            f"catalogs, but got npatch = {npatches}. If these catalogs were " + \
            "loaded from saved treecorr catalogs, check that they were " + \
            "written by the same run."
        )

    print(f"All catalogs share {list(npatches.values())[0]} patches\n")


def get_dust(fg, fgr, bg, bgr, names, correl_config):
    '''
    Run correlations, save to file
    '''

    print('Correlating fg x bg...\n')
    DK = treecorr.NKCorrelation(**correl_config['treecorr_params'])
    DK.process(fg.treecorrCatalog, bg.treecorrCatalog)
    DK.write(names.dk_outfile)

    print('Correlating fg_rand x bg...\n')
    FR = treecorr.NKCorrelation(**correl_config['treecorr_params'])
    FR.process(fgr.treecorrCatalog, bg.treecorrCatalog)
    FR.write(names.fr_outfile)

    print('Calculating compensated signal...\n')
    corr_xi, corr_varxi = DK.calculateXi(rk=FR)
    DK.write(rk=FR, file_name=names.ck_outfile)

    print('Correlating fg x bg_rand...\n')
    RK = treecorr.NKCorrelation(**correl_config['treecorr_params'])
    RK.process(fg.treecorrCatalog, bgr.treecorrCatalog)
    RK.write(names.dr_outfile)

    print('Correlating fg_rand x bg_rand...\n')
    RR = treecorr.NKCorrelation(**correl_config['treecorr_params'])
    RR.process(fgr.treecorrCatalog, bgr.treecorrCatalog)
    RR.write(names.rr_outfile)

    print('Calculating fg/fgr/bg/bgr covariance...\n')
    var_method = correl_config['treecorr_params']['var_method']
    jointcov = treecorr.estimate_multi_cov([DK, RK, RR], var_method)
    np.savetxt(names.cov_output, jointcov)


def main(args):

    # This is the Menard redshift
    z_theory = 0.36

    # Read in configuration file
    config_file = args.config
    correl_config = utils.read_yaml(config_file)

    # Match TreeCorr's OpenMP pool to the cores we actually hold
    set_treecorr_threads(correl_config)

    # Create output directory if it doesn't exist
    if not os.path.isdir(correl_config['output_path']):
        os.makedirs(correl_config['output_path'])

    # Load background catalog & calculate dust reddening
    # This sets patch centers
    bg = Correlator(correl_config, ctype='background_catalog')
    bg.load()
    # If we have read in a TreeCorr catalog from file, these do nothing
    bg.do_reddening()
    bg.write_treecorr_cat_to_file()

    # Grab patch centers once; for a bg read from file they are derived from
    # its patch column, which means reading the whole catalog
    patch_centers = bg.treecorrCatalog.patch_centers

    # Load background random catalog & calculate dust reddening
    # Include background patch_centers for covariance calculations
    bgr = Correlator(correl_config, ctype='background_randoms')
    bgr.load(treecorr_patch_centers=patch_centers)
    # Again, if we have read in a TreeCorr catalog from file, these do nothing
    bgr.do_reddening()
    bgr.write_treecorr_cat_to_file()

    # Load foreground catalog
    fg = Correlator(correl_config, ctype='foreground_catalog')
    fg.load(treecorr_patch_centers=patch_centers)
    mean_fg_z = get_mean_fg_redshift(fg, correl_config)

    # Load foreground random catalog
    fgr = Correlator(correl_config, ctype='foreground_randoms')
    fgr.load(treecorr_patch_centers=patch_centers)

    # Guard against catalogs disagreeing about patches
    check_patch_consistency(correl_config, fg=fg, fgr=fgr, bg=bg, bgr=bgr)

    # Make names
    names = make_names(correl_config)

    if args.saveonly != True: 
        # Do calculation
        get_dust(
            fg=fg, fgr=fgr, bg=bg, bgr=bgr,
            names=names, correl_config=correl_config
        )

        # Make pretty output plots
        print('Plotting output figure...\n')
        plot = DustPlotter(
            dk_file = names.dk_outfile,
            dr_file = names.dr_outfile,
            fr_file = names.fr_outfile,
            rr_file = names.rr_outfile,
            ck_file = names.ck_outfile,
            z_fg = mean_fg_z,
            z_theory = z_theory
        )
        plot.plot_res(outplotn=names.fig_output, kpc=correl_config['use_kpc'])

    else:
        print("TreeCorr catalogs saved to file")
        print("Program complete, exiting.")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Runner script for dust calculation operations."
        )
    parser.add_argument(
        "--config","-c", type=str, required=True,
        help="Path to the configuration file."
        )
    parser.add_argument(
        "--saveonly", action='store_true',
        help="Save TreeCorr catalog outputs only, then exit [default: False]"
        )
    args = parser.parse_args()
    main(args)
