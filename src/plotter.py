import numpy as np
from matplotlib import rc,rcParams
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import healpy as hp
import pdb
import fitsio
from astropy.cosmology import Planck18 as cosmo
import astropy.coordinates as coord
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table


class RCParamsMixin:
    ##
    ## Make plot settings pretty
    ##

    def set_rc_params(self, fontsize=None):

        print("Setting Matplotlib RC parameters...")

        if fontsize is None:
            fontsize=16
        else:
            fontsize=int(fontsize)

        rc('font',**{'family':'serif'})
        rc('text', usetex=True)

        #plt.rcParams.update({'figure.facecolor':'w'})
        plt.rcParams.update({'axes.linewidth': 1.3})
        plt.rcParams.update({'xtick.labelsize': fontsize})
        plt.rcParams.update({'ytick.labelsize': fontsize})
        plt.rcParams.update({'xtick.major.size': 8})
        plt.rcParams.update({'xtick.major.width': 1.3})
        plt.rcParams.update({'xtick.minor.visible': True})
        plt.rcParams.update({'xtick.minor.width': 1.})
        plt.rcParams.update({'xtick.minor.size': 6})
        plt.rcParams.update({'xtick.direction': 'out'})
        plt.rcParams.update({'ytick.major.width': 1.3})
        plt.rcParams.update({'ytick.major.size': 8})
        plt.rcParams.update({'ytick.minor.visible': True})
        plt.rcParams.update({'ytick.minor.width': 1.})
        plt.rcParams.update({'ytick.minor.size':6})
        plt.rcParams.update({'ytick.direction':'out'})
        plt.rcParams.update({'axes.labelsize': fontsize})
        plt.rcParams.update({'axes.titlesize': fontsize})
        plt.rcParams.update({'legend.fontsize': int(fontsize-2)})

        return


class DustPlotter(RCParamsMixin):

    def __init__(self, dk_file=None, dr_file=None, fr_file=None,
                    rr_file=None, ck_file=None, z_fg=None, z_theory=None):
        '''
        cat_files: str
            Filename for correlation functions.
        z_theory: float
            Average value of Menard 2010 galaxies
        z_bg: float
            Mean redshift of background galaxies
        '''

        self.dk = None
        self.dr = None
        self.fr = None
        self.rr = None
        self.ck = None  # Compensated signal
        self.z_fg = z_fg
        self.z_theory = z_theory

        self._load_catalogs(dk_file, dr_file, fr_file, rr_file, ck_file)

    def _load_catalogs(self, dk_file, dr_file, fr_file, rr_file, ck_file):
        '''
        dk_file : the signal (foreground x background galaxies correlation)
        dr_file : foreground x random background
        fr_file : random foreground x background
        rr_file : random foreground x random background
        ck_file : compensated foreground x background signal
        '''
        
        self.dk = Table.read(dk_file, format='ascii', header_start=1)
        self.dr = Table.read(dr_file, format='ascii', header_start=1)
        self.fr = Table.read(fr_file, format='ascii', header_start=1)
        self.rr = Table.read(rr_file, format='ascii', header_start=1)
        if ck_file != None:
            self.ck = Table.read(ck_file, format='ascii', header_start=1)

    def _plot_res(self, outplotn, kpc):
        '''
        Plot up the correlation results, compare to the published result.
        dk : fg x bg correlation
        dr : fg x random
        fr : random fg x bg
        rr : fg random x bg random
        dk_corr: corrected correlation function
        kpc : if true, plot x-axis in kpc units not arcmin (default false)
        outplotn : name for file to plot
        '''
        
        dk = self.dk
        dr = self.dr
        fr = self.fr
        rr = self.rr
        dk_corr = self.ck

        self.set_rc_params()
 
        # Scale Menard theory relationship to present-day
        # Proper kpc yields answer closest to the "book answer"
        theory_kpc = cosmo.kpc_proper_per_arcmin(self.z_theory)
        fg_gal_kpc = cosmo.kpc_proper_per_arcmin(self.z_fg)

        if kpc == True:
            # Converts us from arcminutes to kpc
            scl = fg_gal_kpc.value
            #theory_scl = 304*cosmo.h -- yields their answer
            theory_scl = theory_kpc.value
            xlabel_unit = 'kpc'
        else:
            # Keep plot in arcminutes
            scl = 1.0
            #theory_scl  = (fg_gal_kpc/theory_kpc).value
            theory_scl = theory_kpc.value
            xlabel_unit = 'arcmin'

        # Not sure this is right, trying to go from arcmin at z=0.36 to
        # equivalent arcmin at z=0.11
        # r is defined in arcminutes, so keep it in arcminutes.
        # The relationship does seem to be defined in "inverse arcminutes"/"inverse kpc"
        # so scales need to be 1/kpc_per_arcmin
        # SO, 1/arcmin *(1/(kpc_per_arcmin)) -- arcmin/kpc.
        theory_r = np.logspace(-2,5,10)
        av = 2.5e-3 * (theory_r/theory_scl)**(-0.86)

        fig, ax = plt.subplots(figsize=(10,7), tight_layout=True)

        ax.plot(theory_r, av, color='tab:red',
                label=f'Menard (2010) scaled to z={self.z_fg:.3f}')
        ax.errorbar(dk['meanr']*scl, dk['kappa'], yerr=dk['sigma'],
                    fmt='-.o', capsize=5, color='tab:orange',
                    label='raw signal')
        ax.errorbar(dk['meanr']*scl, dk['kappa']-fr['kappa'],
                    yerr=dk['sigma'], fmt='-o', capsize=5,
                    color='tab:green', label='signal - fg_random')
        ax.errorbar(dk['meanr']*scl,
                    dk['kappa']-fr['kappa']-dr['kappa']+rr['kappa'],
                    yerr=dk['sigma'], fmt='-o', capsize=5, color='tab:blue',
                    label='signal - all_randoms')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.05*scl, 200*scl)
        ax.set_ylim(1E-5, 1)
        ax.set_xlabel(f'Impact parameter ({xlabel_unit})', fontsize=16)
        ax.set_ylabel(r'$A_{\rm V}$ (mag)', fontsize=16)
        ax.set_title('SCOS x redMaGiC', fontsize=16)
        ax.legend(fontsize=14)

        fig.savefig(outplotn)
        
        try:
            fig, ax = plt.subplots(figsize=(10,7), tight_layout=True)
            
            ax.plot(theory_r, av, color='tab:red',
                    label=f'Menard (2010) scaled to z={self.z_fg:.3f}')
            ax.errorbar(dk_corr['meanr']*scl, dk_corr['kappa'],
                        yerr=dk_corr['sigma'], fmt='--o', capsize=5, 
                        color='tab:orange', label='compensated signal')
            ax.errorbar(dk_corr['meanr']*scl,
                        dk_corr['kappa']-dr['kappa']+rr['kappa'],
                        yerr=dk_corr['sigma'],
                        fmt='-o', capsize=5, color='tab:blue', 
                        label='compensated signal - all randoms')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(0.05*scl, 200*scl)
            ax.set_ylim(1E-5, 1)
            ax.set_xlabel(f'Impact parameter ({xlabel_unit})', fontsize=16)
            ax.set_ylabel(r'$A_{\rm V}$ (mag)', fontsize=16)
            ax.set_title('SCOS x redMaGiC', fontsize=16)
            ax.legend(fontsize=14)
            
            # Assuming we have gotten thus far, let's make a new file name
            basename = os.path.basename(outplotn)
            dirname  = os.path.dirname(outplotn)
            comp_fig_name = 'compensated_' + basename
            print(f"compensated figure name is {os.path.join(dirname, comp_fig_name)}")
            fig.savefig(os.path.join(dirname, comp_fig_name))
            
        except:
            raise "problem with correlplot"
                    
    def plot_res(self, outplotn='dust_correlation_plot.png',
                    kpc=False, subsample=False):
        '''
        outplotn : save plot file to
        kpc : if true, plot x-axis in kpc units (default)
        subsample : use special-purpose plotting software
        '''
        
        if subsample == True:
            raise AttributeError('subsample plot is deprecated')
            #self._plot_corr_res(outplotn, kpc)
        else:
            print(f'saving correlation plot figure to {outplotn}')
            self._plot_res(outplotn, kpc)
        

class OverlapPlotter(RCParamsMixin):

    def __init__(self, cat1_name=None, cat2_name=None,
                    outname='catalog_overlap.pdf', outdir='./'):
        '''
        Make nice catalog RA/Dec overlap plot in both Cartesian and
        Aitoff projection

        Inputs
            cat1: path to first catalog
            cat2: path to second catalog
            outname: name to save figure to
            outdir: where should file be saved

        TO DO: add smarter error handling for the cat1/cat2 objects
        '''

        self.cat1_name = cat1_name
        self.cat2_name = cat2_name
        self.outname = outname
        self.outdir = outdir


        if type(cat1_name)==str:
            self.cat1 = Table.read(cat1_name, memmap=True)
        else:
            self.cat1 = cat1_name

        if type(cat2_name)==str:
            self.cat2 = Table.read(cat2_name, memmap=True)
        else:
            self.cat2 = cat2_name


    def make_plot(self, outname=None, projection=None,
                  ra_tag1=None, dec_tag1=None, coordframe1=None, label1=None,
                  ra_tag2=None, dec_tag2=None, coordframe2=None, label2=None):
        '''
        TO DO: need to fix catalog labels, also find a smarter way
        to deal with coordinate systems of the two.

        Note: default "projection" is rectilinear; other options include
        "mollweide", "aitoff", and "hammer" (time).
        '''

        self.set_rc_params(fontsize=16)

        if outname is None:
            outname = self.outname

        if (os.path.dirname(outname)==''):
            outname=os.path.join(self.outdir, outname)

        # Create SkyCoord object to hold RA, Dec of catalogs
        cat1 = self.cat1; cat2 = self.cat2

        # Labels
        if label1 == None:
            label1 = os.path.basename(self.cat1_name)
        if (cat2 is not None) & (label2 == None):
            label2 = os.path.basename(self.cat2_name)

        if (ra_tag1 is not None) & (coordframe1 is not None):
            gal1 = SkyCoord(cat1[ra_tag1], cat1[dec_tag1],
                            frame=coordframe1, unit=u.deg)
            sky1 = gal1.icrs
        else:
            try:
                gal1 = SkyCoord(cat1['l'], cat1['b'],
                                frame='galactic', unit=u.deg)
                sky1 = gal1.icrs
            except KeyError:
                sky1 = SkyCoord(cat1['ra'], cat1['dec'],
                                frame='icrs', unit=u.deg)

        if (ra_tag2 is not None) & (dec_tag2 is not None) \
            & (coordframe2 is not None):
            gal2 = SkyCoord(cat2[ra_tag2], cat2[dec_tag2],
                            frame=coordframe2, unit=u.deg)
            sky2 = gal2.icrs
        else:
            try:
                gal2 = SkyCoord(cat2['l'], cat2['b'],
                                frame='galactic', unit=u.deg)
                sky2 = gal2.icrs
            except KeyError:
                sky2 = SkyCoord(cat2['ra'], cat2['dec'],
                                frame='icrs', unit=u.deg)
            except TypeError:
                sky2 = None

        # Create a plot instance (can also use axes class)
        # Note: aitoff projection apparently avoids mollweide's extreme edge distortions
        if projection in ['mollweide', 'aitoff']:
            figsize=(11.5, 6)
        else:
            figsize=(10, 7)

        fig, ax = plt.subplots(1,1, figsize=figsize, tight_layout=True, \
                        subplot_kw=dict(projection=projection))
        ax.grid(True)
        ax.set_xlabel('RA'); ax.set_ylabel('Dec')

        # Plot the points - it takes a long time for them all to show up!
        ax.plot(sky1.ra.wrap_at('180d').radian, sky1.dec.radian, '.',
                label=label1, color='xkcd:marine', markersize=0.025)
        if (sky2 is not None):
            ax.plot(sky2.ra.wrap_at('180d').radian, sky2.dec.radian, '.',
                    label=label2, color='xkcd:neon red', markersize=0.025)

        ax.legend(markerscale=400, loc='upper right')
        fig.tight_layout()

        fig.savefig(outname)
        fig.savefig(outname.replace('pdf', 'png'))


    def make_Av_map(self, outname=None, projection=None, label1=None):
        '''
        The input should be a treecorrcat;
        Note: default "projection" is rectilinear; other options include
        "mollweide", "aitoff", and "hammer" (time).
        '''

        self.set_rc_params(fontsize=16)

        # Sanity check
        #if ['ra', 'dec', 'k'] not in self.cat1.__dict__.keys():
        #    print("Supplied catalog is missing one of: 'ra', 'dec', 'k'")

        if outname is None:
            outname = self.outname

        if (os.path.dirname(outname)==''):
            outname=os.path.join(self.outdir, outname)

        # Labels
        if label1 == None:
            label1 = os.path.basename(self.cat1_name)

        # Create SkyCoord object to hold RA, Dec of catalogs
        sky1 = SkyCoord(self.cat1['ra'], self.cat1['dec'],
                        frame='icrs', unit=u.deg)

        # Create a plot instance (can also use axes class). Note that aitoff
        # projection apparently avoids mollweide's extreme edge distortions
        if projection in ['mollweide', 'aitoff']:
            figsize=(11.5, 6)
        else:
            figsize=(10, 7)

        fig, ax = plt.subplots(1,1, figsize=figsize, tight_layout=True, \
                        subplot_kw=dict(projection=projection))
        ax.grid(True)
        ax.set_xlabel('RA'); ax.set_ylabel('Dec')

        # Plot the points - it takes a long time for them all to show up!
        ax.scatter(sky1.ra.wrap_at('180d').radian, sky1.dec.radian, ',',
                    label=label1, c=self.cat1['k'].data)
        ax.legend(markerscale=400, loc='upper right')
        fig.tight_layout()

        fig.savefig(outname)
        fig.savefig(outname.replace('pdf', 'png'))
