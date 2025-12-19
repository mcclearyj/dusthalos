###
### Module to calculate source density around lenses using TreeCorr
### Technically, this isn't a source density calculation in terms of 
### number counts but rather a cross-correlation between lens positions
### and source positions. Still useful for quantifying any change in 
### source/bg density as a function of angular separation from lenses/fg.
###

import treecorr
import numpy as np
import os
import pdb  # For debugging

### Local imports
import src.plotter as plotter
import src.utils as utils
from src.correlator import Correlator

def make_nn(var_method='shot'):      
    return treecorr.NNCorrelation(min_sep=theta_edges[0], max_sep=theta_edges[-1],
                                  nbins=len(theta_edges)-1, sep_units='arcmin', 
                                  bin_slop=0, var_method=var_method)

def calculate_xi(fg, fgr, bg, bgr, var_method=None):
    """
    Run the correlation calculation with the specified variance method. 
    Input: fg, bg, fgr, bgr are instances of dusthalos.Catalog(), which include 
    TreeCorr Catalog instances as attributes.
    Returns the xi and varxi arrays.
    """
    if type(fg) is not treecorr.Catalog or type(bg) is not treecorr.Catalog or \
        type(fgr) is not treecorr.Catalog or type(bgr) is not treecorr.Catalog:
        raise TypeError("Input catalogs must be TreeCorr Catalog objects")

    if var_method is None: 
        print("No covariance method specified, defaulting to shot noise only")
        var_method = 'shot'

    dd = _calculate_xi(fg, bg, fgr, bgr, var_method=var_method)
    return dd

def _calculate_xi(fg, bg, fgr, bgr, var_method='shot'):
    """
    This is essentially a Landy-Szalay estimator clustering calculation.
    Note that input catalogs are assumed to be TreeCorr Catalog objects.
    """
    # Intitialize TreeCorr NNCorrelation objects for each of the four catalogs
    dd = make_nn(var_method=var_method); dr = make_nn(var_method=var_method) 
    rd = make_nn(var_method=var_method); rr = make_nn(var_method=var_method)

    print("Processing TreeCorr correlations...")
    dd.process(fg,  bg)    # D_l D_s
    dr.process(fg,  bgr)   # D_l R_s

    # BREAK FOR SANITY CHECK! This is an approximation for the boot factor (B)
    # Normalize by catalog sizes the way TreeCorr does internally:
    # Pair counts scale ~ N_l * N_s, so correct the ratio by N_Rs / N_Ds
    # In this case, N_Rs / N_Ds == len(bgr.ra) / len(bg.ra)
    B = (dd.npairs / dr.npairs) * (len(bgr.ra) / len(bg.ra))
    print("Quick-check boost factor B= ", B)

    # Back to business. Process the other two combinations
    rd.process(fgr, bg)    # R_l D_s
    rr.process(fgr, bgr)   # R_l R_s

    print("Proceeding with full clustering/xi calculation...")
    xi, varxi = dd.calculateXi(rr=rr, dr=dr, rd=rd)
    B = 1.0 + xi

    # Printing out results here since TreeCorr has some funny behavior 
    # when things get written to file. 
    print("xi: ", xi)
    print("Boost factor B: ", B)
    print("varxi: ", varxi)
    return dd

def load_catalogs(correl_config):
    """
    Use Correlator class to load the necessary catalogs, including 
    conversion of coordinates to ICRS format, and create instances of 
    TreeCorr catalogs for each. 
    """
    # Load background catalog and set patch centers defined in correl_config
    bg = Correlator(correl_config, ctype='background_catalog')
    bg.load()

    # Load background random catalog
    # Include background patch_centers for covariance calculations
    bgr = Correlator(correl_config, ctype='background_randoms')
    bgr.load(treecorr_patch_centers=bg.treecorrCatalog.patch_centers)

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

    return fg, bg, fgr, bgr

def run_all(correl_config, theta_edges=None):
    """
    Run the correlation calculation.
    NOTE: this function overrides the correl_config angular separation
    bin values, since this source density calculation is special-purpose.
    """

    # Define theta_edges (angular separations in arcminutes) for TreeCorr
    if theta_edges is None:
        theta_edges = np.array([0.5, 1, 2, 5, 10, 20, 50, 100, 200])
    print("Using default theta edges:", theta_edges)

    # Load catalogs and return dusthalos.Catalog() instances, 
    # which include TreeCorr Catalog instances as attributes
    print("Loading catalogs...")
    fg, bg, fgr, bgr = load_catalogs(correl_config)

    # Calculate xi and varxi aka clustering
    print("Calculating xi and varxi...")
    dd = calculate_xi(
        fg.treecorrCatalog, bg.treecorrCatalog, 
        fgr.treecorrCatalog, bgr.treecorrCatalog, 
        var_method=correl_config['treecorr_params']['var_method']
    )

    # Write output to file
    print("Writing correlation result to file...")
    dd.write(os.path.join(
        correl_config['output_dir'], correl_config['output_basename']+'clustering.txt'
    ))

    print("Clustering calculation complete.")


