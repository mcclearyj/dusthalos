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

    # Create output directory if it doesn't exist
    if not os.path.isdir(correl_config['output_path']):
        os.makedirs(correl_config['output_path'])

    # Load background catalog & calculate dust reddening
    # This sets patch centers
    bg = Correlator(correl_config, ctype='background_catalog')
    bg.load()
    bg.do_reddening()
    bg.write_to_file()

    # Load background random catalog & calculate dust reddening
    # Include background patch_centers for covariance calculations
    bgr = Correlator(correl_config, ctype='background_randoms')
    bgr.load(treecorr_patch_centers=bg.treecorrCatalog.patch_centers)
    bgr.do_reddening()

    # Load foreground catalog
    fg = Correlator(correl_config, ctype='foreground_catalog')
    fg.load(treecorr_patch_centers=bg.treecorrCatalog.patch_centers)
    try:
        mean_fg_z = np.median(fg.Catalog.data[fg.cat_config['z_key']])
    except KeyError:
        print("No redshift column found, setting mean fg redshift to 0")
        print("Forcing arcminute plot scaling")
        mean_fg_z = 0; correl_config['use_kpc'] = False

    # Load foreground random catalog
    fgr = Correlator(correl_config, ctype='foreground_randoms')
    fgr.load(treecorr_patch_centers=bg.treecorrCatalog.patch_centers)

    # Make names
    names = make_names(correl_config)

    # Do calculation
    get_dust(fg=fg, fgr=fgr, bg=bg, bgr=bgr,
                names=names, correl_config=correl_config)

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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Runner script for Catalog operations."
    )
    parser.add_argument(
        "-config","-c", type=str,
        help="Path to the configuration file.", required=True
    )
    args = parser.parse_args()
    main(args)
