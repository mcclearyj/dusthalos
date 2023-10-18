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
                    rr_file=None, z_fg=None, z_theory=None):
        '''
        cat_files: str
            Filename for correlation functions.
        z_theory: float
            Average value of Menard 2010 galaxies
        z_bg: float
            Mean redshift of background galaxies
        '''

        self.dk = dk_file
        self.dr = dr_file
        self.fr = fr_file
        self.rr = rr_file
        self.z_fg = z_fg
        self.z_theory = z_theory

        self._load_catalogs(dk_file, dr_file, fr_file, rr_file)

        return


    def _load_catalogs(self, dk_file, dr_file, fr_file, rr_file):
        '''
        dk_file : the signal (foreground x background galaxies correlation)
        dr_file : foreground x random background
        fr_file : random foreground x background
        rr_file : random foreground x random background
        '''

        if dk_file is None:
            dk_file = 'dust_correlation_signal.fits'
        if dr_file is None:
            dr_file = 'dust_correlation_bg_randoms.fits'
        if fr_file is None:
            fr_file = 'dust_correlation_fg_randoms.fits'
        if rr_file is None:
            rr_file = 'dust_correlation_fgxbg_randoms.fits'

        self.dk = fitsio.read(dk_file)
        self.dr = fitsio.read(dr_file)
        self.fr = fitsio.read(fr_file)
        self.rr = fitsio.read(rr_file)

        return


    def _plot_res(self, outplotn, kpc):
        '''
        Plot up the correlation results, compare to the published result.
        dk : fg x bg correlation
        dr : fg x random
        fr : random fg x bg
        rr : fg random x bg random
        kpc : if true, plot x-axis in kpc units not arcmin (default false)
        outplotn : name for file to plot
        '''

        dk = self.dk
        dr = self.dr
        fr = self.fr
        rr = self.rr

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
            theory_scl  = (fg_gal_kpc/theory_kpc).value
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

        ax.plot(theory_r,av,label=f'Menard (2010) scaled to z={self.z_fg:.3f}', color='tab:red')


        ax.errorbar(dk['meanr']*scl, dk['kappa'], yerr=dk['sigma'],
                    fmt='-.o',capsize=5, color='tab:orange',
                    label='raw signal')
        ax.errorbar(dk['meanr']*scl, dk['kappa']-fr['kappa'],
                    yerr=dk['sigma'], fmt='-',capsize=5,
                    color='tab:green', label='signal - fg_random')
        ax.errorbar(dk['meanr']*scl,
                    dk['kappa']-fr['kappa']-dr['kappa']+rr['kappa'],
                    yerr=dk['sigma'],fmt='-o',capsize=5, color='tab:blue',
                    label='signal - all_randoms')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.05*scl, 200*scl)
        ax.set_ylim(1E-5, 1)
        ax.set_xlabel(f'Impact parameter ({xlabel_unit})', fontsize=16)
        ax.set_ylabel(r'$A_v$ (mag)', fontsize=16)
        ax.set_title('SCOS x redMaGiC', fontsize=16)
        ax.legend(fontsize=14)

        fig.savefig(outplotn)


        return


    def _plot_res_subsample(self):
        '''
        Special purpose plotting code for the report, saved b/c fancy error bar

        dk : fg x bg correlation
        dr : fg x random
        fr : random fg x bg
        rr : fg random x bg random
        outplotn : name for file to plot
        '''

        dk = self.dk
        dr = self.dr
        fr = self.fr
        rr = self.rr

        low_dk = dk[0:3]
        hi_dk = dk[3:]
        mean_lo_dk = np.mean(low_dk['kappa'])
        mean_lo_r = np.mean(low_dk['meanr'])
        mean_lo_var = low_dk['sigma'][-1]

        kmr = dk['kappa']-fr['kappa']-dr['kappa']+rr['kappa']
        low_kmr = kmr[0:3]
        mean_lo_kmr = np.mean(low_kmr)


        fake_r = [mean_lo_r]; fake_r.extend(dk['meanr'][3:])
        fake_dk = [mean_lo_dk+mean_lo_var]; fake_dk.extend(dk['kappa'][3:5]+dk['sigma'][3:5]); fake_dk.extend(dk['kappa'][5:])
        #fake_kmr = [mean_lo_kmr+mean_lo_var]; fake_kmr.extend(kmr['kappa'][3:5]+kmr['sigma'][3:5]); fake_kmr.extend(dk['kappa'][5:])
        fake_kmr = [mean_lo_kmr+mean_lo_var]; fake_kmr.extend(kmr[3:6]+dk['sigma'][3:6]);fake_kmr.extend(kmr[6:])
        fake_sig = [low_dk['sigma'][-1]]; fake_sig.extend(dk['sigma'][3:])
        fake_dk1 = [mean_lo_dk+mean_lo_var]; fake_dk1.extend(dk['kappa'][3:])

        fake_r = np.array(fake_r)
        dk = np.array(dk)
        scl = 2.017*60 / 1000
        #dk['meanr']*=scl
        #fake_r*=scl

        fig = plt.figure(figsize=(10,7), tight_layout=True)
        ax=fig.add_subplot(111)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.plot(fake_r[0:6]*scl, fake_dk1[0:6], ':o', lw=1.2, color='tab:orange',  label=None)
        ax.plot(fake_r[0:6]*scl, fake_kmr[0:6], ':o', lw=1.2, color='tab:blue',  label=None)
        ax.plot(fake_r[5:]*scl, fake_dk1[5:], '-o', lw=1., color='tab:orange',  label='raw signal')
        ax.plot(fake_r[5:]*scl, fake_kmr[5:], '-o', lw=1., color='tab:blue',  label='signal -  randoms')


        err0, caps, bars = ax.errorbar(mean_lo_r*scl,mean_lo_dk+mean_lo_var, yerr=mean_lo_var, color='tab:orange',  capsize=5, uplims=True,label=None)
        caps[1].set_marker(',')
        err3, caps, bars = ax.errorbar(dk['meanr'][3]*scl, dk['kappa'][3], yerr=dk['sigma'][3]*0.65, capsize=5,color='tab:orange', uplims=True,label=None)
        caps[1].set_marker(',')

        err4 = ax.errorbar(dk['meanr'][4]*scl, dk['kappa'][4], yerr=dk['sigma'][4], capsize=5,color='tab:orange',label=None)
        err5,caps,bar = ax.errorbar(dk['meanr'][5]*scl, dk['kappa'][5], yerr=dk['sigma'][5], lolims=True, capsize=5,color='tab:orange', label=None)
        caps[1].set_marker(',')

        err6 = ax.errorbar(dk['meanr'][6:]*scl, dk['kappa'][6:], yerr=dk['sigma'][6:], fmt='o', capsize=5,color='tab:orange', label=None)

        err0,caps,bars = ax.errorbar(mean_lo_r*scl,mean_lo_kmr+mean_lo_var, yerr=mean_lo_var, color='tab:blue',  capsize=5, uplims=True,label=None)
        caps[1].set_marker(',')

        err3,caps,bars = ax.errorbar(dk['meanr'][3]*scl, kmr[3]+dk['sigma'][3], yerr=dk['sigma'][3], capsize=5,color='tab:blue',uplims=True, label=None)
        caps[1].set_marker(',')

        err4,caps,bars = ax.errorbar(dk['meanr'][4]*scl, kmr[4]+dk['sigma'][4], yerr=dk['sigma'][4]*0.9, capsize=5,color='tab:blue',uplims=True, label=None); caps[1].set_marker(',')
        err5,caps,bars = ax.errorbar(dk['meanr'][5]*scl, kmr[5]+dk['sigma'][5], yerr=dk['sigma'][5]*0.6, capsize=5,color='tab:blue',uplims=True, label=None);caps[1].set_marker(',')
        err6=ax.errorbar(dk['meanr'][6:]*scl, kmr[6:], yerr=dk['sigma'][6:], capsize=5,color='tab:blue', label=None)

        #ax.plot(dk['meanr'],dk['kappa'],'o' , color='tab:orange', label='raw signal')
        #ax.plot(dk['meanr'],kmr,'o', color='tab:blue', label='signal - randoms')

        #ax.errorbar
        #ax.errorbar(dk['meanr'][4:], dk['kappa'][4:], yerr=dk['sigma'][4:], fmt='o',linestyle=':', capsize=5,color='tab:orange', label='raw signal')
        #ax.errorbar(dk['meanr'][4:], kmr[4:], yerr=dk['sigma'][4:],fmt='o', linestyle='-', color='tab:blue',capsize=5, label='signal - randoms')

        ax.set_ylim(1e-4,0.4)
        ax.set_xlim(0.026,25)
        r = np.logspace(-2,2.7,10)*scl
        av = 2.4e-3 * (r/2.227)**(-0.84)
        avplot = ax.plot(r,av,label='Menard+ 2010 (scaled)',color='tab:red')
        ax.set_xlabel('Impact parameter (kpc)',fontsize=16)
        ax.set_ylabel(r'$A_v$ (mag)',fontsize=16)
        ax.set_title(r'SCOS $\times$ redMaGiC', fontsize=16)
        ax.legend(fontsize=14)

        fig.tight_layout()
        fig.savefig(outplotn)

        return


    def plot_res(self, outplotn='dust_correlation_plot.png',
                    kpc=False, subsample=False):
        '''
        outplotn : save plot file to
        kpc : if true, plot x-axis in kpc units (default)
        subsample : use special-purpose plotting software
        '''

        if subsample == True:
            self._plot_res_subsample()

        else:
            print(f'saving correlation plot figure to {outplotn}')
            self._plot_res(outplotn, kpc)

        return


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
                    markersize=0.02, label=label1, color='xkcd:bluegrey')
        if (sky2 is not None):
            ax.plot(sky2.ra.wrap_at('180d').radian, sky2.dec.radian, '.',
                        markersize=0.02, label=label2, color='xkcd:neon red')
        ax.legend(markerscale=400, loc='upper right')
        fig.tight_layout()

        fig.savefig(outname)
        fig.savefig(outname.replace('pdf', 'png'))
