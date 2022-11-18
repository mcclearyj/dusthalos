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

def rc_params(self):

    rc('font',**{'family':'serif'})
    rc('text', usetex=True)

    rcParams['axes.linewidth'] = 1.3
    rcParams['xtick.labelsize'] = 16
    rcParams['ytick.labelsize'] = 16

    rcParams['xtick.minor.visible'] = True
    rcParams['xtick.major.width'] = 1.2
    rcParams['xtick.minor.width'] = 1.2

    rcParams['xtick.direction'] = 'out'
    rcParams['ytick.minor.visible'] = True
    rcParams['ytick.major.width'] = 1.2
    rcParams['ytick.direction'] = 'out'
    rcParams['ytick.minor.width'] = 1.2

return


class DustPlotter():

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

        self.dk = None
        self.dr = None
        self.fr = None
        self.rr = None
        self.z_fg = None
        self.z_theory = None

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


    def _plot_res(self, outplotn):
        '''
        Plot up the correlation results, compare to the
        published result.

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

        # Convert from arcmin to kpc
        theory_scl = cosmo.kpc_proper_per_arcmin(self.z_theory)
        scl = cosmo.kpc_proper_per_arcmin(self.z_fg)

        r = np.logspace(-2,2.7,10)*scl
        # Converts us from arcminutes to kpc
        av = 2.5e-3 * (r)**(-0.86) # Ideally, will do a shaded plot with error bar of relationship

        fig, ax = plt.subplots(figsize=(10,7), tight_layout=True)
        ax.set_xscale('log')
        ax.set_yscale('log')



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
        ax.set_xlabel('Impact parameter (Mpc)',fontsize=16)
        ax.set_ylabel(r'$A_v$ (mag)',fontsize=16)
        ax.set_title(r'SCOS $\times$ redMaGiC', fontsize=16)
        ax.legend(fontsize=14)

        fig.tight_layout()
        fig.savefig(outplotn)

        return

        def run(self, z_fg=None, z_th=0.36, )
