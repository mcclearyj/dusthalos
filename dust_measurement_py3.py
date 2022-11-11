import numpy as np
from numpy.linalg import norm
from matplotlib import rc,rcParams
rc('font',**{'family':'serif'})
rc('text', usetex=True)

#matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ion()
from astropy.io import fits
import treecorr
import os
from astropy.table import Table, vstack, hstack
import healpy as hp
import pdb
import fitsio


def hpRaDecToHEALPixel(ra, dec, nside=  4096, nest= False):
    phi = ra * np.pi / 180.0
    theta = (90.0 - dec) * np.pi / 180.0
    hpInd = hp.ang2pix(nside, theta, phi, nest= nest)
    return hpInd


def get_fg_catalog(fg_file, vb=False):
    '''
    Load foreground catalog, return a treecorr.Catalog() object
    '''

    if vb is True:
        print(f'loading file {fg_file}...')
    data=Table.read(fg_file,format='fits')
    ra_cat =  data['ra']
    dec_cat = data['dec']

    if vb is True:
        print(f'making treecorr.Catalog for {fg_file}')
    fgCat = treecorr.Catalog(ra=ra_cat, dec=dec_cat, \
            ra_units='deg',dec_units='deg')
    print("Length of catalog = %i" % len(ra_cat))

    return fgCat


def get_fg_randoms(maskfile, nrand=1e6, nside=4096, partial=True, nest=False, vb=False):
    # Make randoms on the sphere.
    if vb == True:
        print("loading hmap, making randoms")
    hmap, hd = hp.read_map(maskfile, partial=partial, nest=nest, h=True)
    # But how many? Try to get approximately nrand, if possible.
    fcover = np.sum(hmap > 0)*1./hmap.size
    ndraw = np.ceil(nrand/fcover*1.2).astype(int)

    ran1, ran2 = np.random.random(2*ndraw).reshape(2, -1)
    ra  = 2*np.pi * (ran1 - 0.5) * 180/np.pi
    dec = np.arcsin(2.*(ran2-0.5)) * 180/np.pi

    hpInd = hpRaDecToHEALPixel(ra,dec,nside=nside,nest=nest)
    keep  = hmap != hp.UNSEEN

    use = np.random.rand(ra.size) < hmap[hpInd]

    ra  = ra[use]
    dec = dec[use]
    ra = ra[:int(nrand)]
    dec = dec[:int(nrand)]

    # covert ra from [-180,180 )  to [0,360)
    # Need to double-check this RA/Dec conversion
    if vb == True:
        print("making fg random treecorr catalog")
    ra = (ra + 360) % 360
    fgRan = treecorr.Catalog(ra=ra, dec=dec, \
            ra_units='deg',dec_units='deg')

    return fgRan


def est_reddening(catalog, ortho=False, ortho_index = 0):
    '''
    Make the colors.
    '''
    
    gmag = catalog['mof_cm_mag_corrected_g']
    rmag = catalog['mof_cm_mag_corrected_r']
    imag = catalog['mof_cm_mag_corrected_i']
    zmag = catalog['mof_cm_mag_corrected_z']

    data = np.vstack([gmag,rmag,imag,zmag])
    covar = np.cov(data)

    # I've forgotten the meaning of these different dust relations.
    if ortho:
        #vdust = np.array([ 4.60, 3.10, 2.18, 1.6 ])
        # Fitzpatrick99
        # vdust = np.array([1.12150099, 0.77164321, 0.57725486, 0.45259124])
        # Calzetti:
        #vdust = np.array([1.13552323, 0.77956032, 0.55082914, 0.39834168])
        # Cardelli, Clayton & Mathis '89
        vdust = np.array([1.12224688, 0.82747095, 0.62680647, 0.47880753])
        dmdp = get_ortho(vdust, index = ortho_index)
        print(f"Dust vector is: {vdust}")
        print(f"ortho dust vector is: {dmdp}")
        print(f"overlap: {np.dot(vdust,dmdp)}")

    else:
        #dmdp = np.array([ 4.60, 3.10, 2.18, 1.6 ])
        # Fitzpatrick99
        # vdust = np.array([1.12150099, 0.77164321, 0.57725486, 0.45259124])
        # Calzetti:
        # vdust = np.array([1.13552323, 0.77956032, 0.55082914, 0.39834168])
        # Cardelli, Clayton & Mathis '89
        dmdp = np.array([1.12224688, 0.82747095, 0.62680647, 0.47880753])

    delta = (np.zeros_like(data).T + dmdp)
    Cinv = np.linalg.inv(covar)
    colors = (data.T - np.average(data,axis=1)).T
    est = np.sum(delta.T*np.dot(Cinv,colors),axis=0)/np.sqrt(np.dot(dmdp,np.dot(Cinv,dmdp)))
    wt = 1./np.sqrt(np.dot(dmdp,np.dot(Cinv,dmdp)))

    return est, wt


def get_bg_catalog(rmz_file, nbins=10, zmin=0.2, vb=False, ortho=False):
    '''
    :rmz_joined: joined catalog with both (MOF) photometry and redshift
    :zmin:       minimum redshift of galaxies to consider
    :nbins:      number of bins for color estimation

    returns treecorr Catalog() object for correlation

    Note: code will break if the mof columns below aren't present
    '''

    if vb == True:
        print(f'Loaded in redMaGic+photometry catalog {rmz_file}')
    try:
        rmz_joined = fitsio.read(rmz_file)
    except:
        rmz_joined = Table.read(rmz_file)

    keep = (rmz_joined['mof_flags']==0) & (rmz_joined['flags_badregions']==0) & \
        (rmz_joined['mof_cm_flux_g']>0) & (rmz_joined['mof_cm_flux_r']>0) & \
        (rmz_joined['mof_cm_flux_i']>0) & (rmz_joined['mof_cm_flux_z']>0) & \
        (rmz_joined['zredmagic'] > zmin)

    # This is a weird digitize thing?
    if vb == True:
        print(f'Beginning background catalog reddening calculation...')

    nbins = nbins
    est = np.zeros(len(rmz_joined))
    est_weight = np.zeros(len(rmz_joined))
    zbins = np.percentile(rmz_joined['zredmagic'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.
    # This feels dangerous: could get messy and ra/dec/z get confused.
    for i in range(nbins):
        these = (rmz_joined['zredmagic'] > zbins[i]) & (rmz_joined['zredmagic'] <= zbins[i+1])
        this_est,this_wt = est_reddening(rmz_joined[these],ortho=ortho)
        est[these] = this_est
        est_weight[these] = this_wt

    if vb == True:
        print(f'Creating background treecorr Catalog...')
    bgCat = treecorr.Catalog(ra=rmz_joined['ra_1'], dec=rmz_joined['dec_1'], \
            k=est, ra_units='deg', dec_units='deg', w=est_weight)

    bgCat.zz = rmz_joined['zredmagic']

    return bgCat


def get_bg_randoms(bg_rand_file, bgCat, nbins=10, zmin=0.2, vb=False):
    '''
    Make random catalog

    NOTE: this is not quite right, as the two different redMaGiC catalogs
    use different random files
    '''
    if vb == True:
        print(f'Reading in background randoms catalog {bg_rand_file}...')
    rm_rand_all = fitsio.read(bg_rand_file)
    rm_rand = rm_rand_all[rm_rand_all['z']>zmin]

    # Do the reddening estimate in redshift slices.
    nbins = nbins
    est = np.zeros(len(rm_rand))
    est_weight = np.zeros(len(rm_rand))
    zbins = np.percentile(rm_rand['z'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.
    # This feels dangerous: could get messy and ra/dec/z get confused.
    for i in range(nbins):
        these = (rm_rand['z'] > zbins[i]) & (rm_rand['z'] <= zbins[i+1])
        these_est = (bgCat.zz > zbins[i]) & (bgCat.zz <=zbins[i+1])
        est[these] = np.random.choice(bgCat.k[these_est],np.sum(these))
        est_weight[these] = np.random.choice(bgCat.w[these_est],np.sum(these))

    if vb == True:
        print('Making bg random treecorr catalog...')
    bgRan = treecorr.Catalog(ra=rm_rand['ra'],dec=rm_rand['dec'],k=est,w=est_weight,\
                                   ra_units='deg',dec_units='deg')

    return bgRan

def plotres(dd_out, dr_out, fr_out=None, rr_out = None, outplotn='fig.png'):
    # Now make a plot.

    dd_outfile = 'dust_correlation_signal.fits'
    dr_outfile = 'dust_correlation_bg_randoms.fits'
    fr_outfile = 'dust_correlation_fg_randoms.fits'
    rr_outfile = 'dust_correlation_fgxbg_randoms.fits'

    dk = fitsio.read(dd_outfile)
    dr = fitsio.read(dr_outfile)
    fr = fitsio.read(fr_outfile)
    rr = fitsio.read(rr_outfile)

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


def main(argv):
    zmin = 0.16
    ortho = False
    vb = True
    nbins = 7
    datapath = '/home/j.mccleary/dust_halos'
    fg_name  = 'SCOS_rmz_match_z_lt_0.15.fits'
    rmz_cat  = 'rmz_y3mof_subset_stacked.fits'
    rm_rand  = 'y3a2_gold2.2.1_redmagic_highdens_randoms.fits'
    rm_mask  = 'y3_gold_2.2.1_RING_joint_redmagic_v0.5.1_wide_maglim_v2.2_mask.fits'

    dd_outfile = 'dust_correlation_signal.fits'
    dr_outfile = 'dust_correlation_bg_randoms.fits'
    fr_outfile = 'dust_correlation_fg_randoms.fits'
    rr_outfile = 'dust_correlation_fgxbg_randoms.fits'
    fig_output = 'dust_correlation_scosRM.png'

    rmz_file = os.path.join(datapath, rmz_cat)
    rmm_file = os.path.join(datapath, rm_mask)
    ran_file = os.path.join(datapath, rm_rand)
    fg_file  = os.path.join(datapath, fg_name)

    print('Getting fg catalog and randoms...\n')
    fgCat = get_fg_catalog(fg_file, vb=vb)
    fgRan = get_fg_randoms(maskfile=rmm_file, partial=True, nest=False)

    print('Getting bg science catalog and randoms...\n')
    bgCat = get_bg_catalog(rmz_file, zmin=zmin, vb=vb, nbins=nbins)
    bgRan = get_bg_randoms(ran_file, bgCat, zmin=zmin, vb=vb, nbins=nbins)

    # Now make the correlation objects. TO DO: make a separate function
    print('Correlating fg x bg...\n')
    DK = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=.6,sep_units='arcmin')
    DK.process(fgCat,bgCat)
    DK.write(dd_outfile)
    print('Correlating fg x bg_rand...\n')
    RK = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=.6,sep_units='arcmin')
    RK.process(fgCat,bgRan)
    RK.write(dr_outfile)
    print('Correlating fg_rand x bg...\n')
    FR = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=.6,sep_units='arcmin')
    FR.process(fgRan,bgCat)
    FR.write(fr_outfile)
    print('Correlating fg_rand x bg_rand...\n')
    RR = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=0.6,sep_units='arcmin')
    RR.process(fgRan,bgRan)
    RR.write(rr_outfile)

    print('Plotting output figure...\n')
    plotres(dd_outfile, dr_outfile, fr_out=fr_outfile, rr_out=rr_outfile, outplotn=fig_output)

    print('Done!')



if __name__ == "__main__":
    import pdb, traceback, sys
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)
