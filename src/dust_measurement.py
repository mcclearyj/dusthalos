import numpy as np
from numpy.linalg import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fitsio
import treecorr
import os
from astropy.table import Table
import pandas as pd
import healpy as hp
import pdb

def get_fg_catalog(datapath,fg_file):
    data=fitsio.read(fg_file)
    try:
        ra_cat =  data['RA']
        dec_cat = data['DEC']
    except:
        try:
            ra_cat =  data['ra']                                                                                                                   
            dec_cat = data['dec']
        except:
            ra_cat = data['_RAJ2000']
            dec_cat = data['_DEJ2000']
       
    catalog = treecorr.Catalog(ra=ra_cat,dec=dec_cat,ra_units='deg',dec_units='deg')
    print( "Length of catalog = %i" % len(ra_cat))
    return catalog

def get_ONbasis(vdust):
    # construct an orthonormal set of basis vectors according
    # to the Gram-Schmidt procedure: define V={v0,v1,v2,v3}
    # that spans the "extinction" basis

    # right now, we have a single vector along which (we think) reddening occurs
    # so define ON basis vectors along which we can project the (real?) vector

    # start by normalizing the input dust vector
    vec = vdust
    
    # Go through the G-S process
    v0 = np.zeros_like(vec); v0[0]= 1.0
    v0prime = v0 - np.dot(v0,vec)*vec
    u0 = v0prime/norm(v0prime)

    v1 = np.zeros_like(vec); v1[1]= 1.0
    v1prime = v1 - np.dot(v1,vec)*vec - np.dot(v1,u0)*u0
    u1 = v1prime/norm(v1prime)

    v2 = np.zeros_like(vec); v2[2]= 1.0
    v2prime = v2 - np.dot(v2,vec)*vec - np.dot(v2,u0)*u0 - np.dot(v2,u1)*u1
    u2 = v2prime/norm(v2prime)

    # Return basis
    return vec,u0,u1,u2
    



def est_reddening(catalog,zeropoint = 30.0, basisvector=None):
    # Make the colors.
    
    gmag = zeropoint - 2.5*np.log10(catalog['mof_flux_g'])
    rmag = zeropoint - 2.5*np.log10(catalog['mof_flux_r'])
    imag = zeropoint - 2.5*np.log10(catalog['mof_flux_i'])
    zmag = zeropoint - 2.5*np.log10(catalog['mof_flux_z'])
    
    data = np.vstack([gmag,rmag,imag,zmag])
    covar = np.cov(data)
    dmdp = basisvector
        
    delta = (np.zeros_like(data).T + dmdp)
    Cinv = np.linalg.inv(covar)
    colors = (data.T - np.average(data,axis=1)).T
    est = np.sum(delta.T*np.dot(Cinv,colors),axis=0)/(np.dot(dmdp,np.dot(Cinv,dmdp)))
    wt = 1./(np.dot(dmdp,np.dot(Cinv,dmdp)))
    
    return est,wt


def get_bg_catalog2(datapath,phot_file,rmz_file,zmin=0.15):

    """ 
    Utility function for reading in the DES background catalog
    Use when the matching catalog has already been made.
    """
    #    joint = os.path.join(datapath,'y1a1_mof_rmz.fits')
    joint = os.path.join(datapath,'redMaGiC_hiz.fits')
    if not os.path.exists(joint):
        bg_phot = Table.read(phot_file,format='fits')
        bg_rmz = Table.read(rmz_file,format='fits')
        print "phot and z files read in"
        # Match these.
        inds = np.in1d(bg_phot['coadd_objects_id'],bg_rmz['ID'])
        cat = bg_phot[inds]
        cat.sort(['coadd_objects_id'])
        bg_rmz.sort(['ID'])
        print "cats matched"
        # Filter these.
        keep = (cat['mof_flags']==0) & (cat['flags_badregion']==0) & \
          (cat['mof_flux_g']>0) & (cat['mof_flux_r']>0) & \
          (cat['mof_flux_i']>0) & (cat['mof_flux_z']>0) & \
          (bg_rmz['ZREDMAGIC']>zmin)
        cat = cat[keep]
        zcat = bg_rmz[keep]
        print "cats filtered"
        # Save to file so we don't have to do this over and over again
        cat.add_columns([zcat['ZREDMAGIC'],zcat['ZREDMAGIC_E']])
        cat.write(joint)
        print "joint cat saved"
    else:
        print("found a joint catalog")
        redMaGiC=fitsio.read(joint,format='fits')
        redMaGiC_filt=redMaGiC[redMaGiC['ZREDMAGIC']>zmin]
        print "background redMaGiC catalog acquired"
    return redMaGiC_filt

def do_reddening_calculation(cat,basisVector):
    print "generating bg catalog for correlation for this vector in ON basis"

    # Do the reddening estimate in redshift slices.
    # Also, loop over vectors in our fake basis

    nbins = 5
    est = np.zeros(cat.size)
    est_weight = np.zeros(cat.size)

    zbins = np.percentile(cat['ZREDMAGIC'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.

    for i in range(nbins):
            these = (cat['ZREDMAGIC'] > zbins[i]) & (cat['ZREDMAGIC'] <= zbins[i+1])
            this_est,this_wt = est_reddening(cat[these],basisvector = basisVector)
            est[these] = this_est
            est_weight[these] = this_wt

    catalog = treecorr.Catalog(ra=cat['ra'],dec=cat['dec'],k=est,\
         ra_units='deg',dec_units='deg',w=est_weight)

    catalog.zz = cat['ZREDMAGIC']
    
    return catalog

def hpRaDecToHEALPixel(ra, dec, nside=  4096, nest= False):
    phi = ra * np.pi / 180.0
    theta = (90.0 - dec) * np.pi / 180.0
    hpInd = hp.ang2pix(nside, theta, phi, nest= nest)
    return hpInd

def get_fg_randoms(nrand = 1e6,maskfile = None,nside=4096,nest=False):
    # Make randoms on the sphere.
    hmap = hp.read_map(maskfile,nest=False)
    # But how many? Try to get approximately nrand, if possible.
    fcover = np.sum(hmap > 0)*1./hmap.size
    ndraw = np.ceil(nrand/fcover*1.2).astype(int)
    
    
    ran1, ran2 = np.random.random(2*ndraw).reshape(2, -1)
    ra = 2*np.pi * (ran1 - 0.5) * 180/np.pi
    dec= np.arcsin(2.*(ran2-0.5)) * 180/np.pi

    
    hpInd = hpRaDecToHEALPixel(ra,dec,nside=nside,nest=nest)
    keep = hmap != hp.UNSEEN
    
    use = np.random.rand(ra.size) < hmap[hpInd]
    ra = ra[use]
    dec= dec[use]
    ra = ra[:np.int(nrand)]
    dec = dec[:np.int(nrand)]

    # covert ra from [-180,180 )  to [0,360)
    ra = (ra + 360) % 360
    rancat = treecorr.Catalog(ra=ra,dec=dec,ra_units='deg',dec_units='deg')
    
    return rancat


def get_bg_randoms(bg_file,Cat,zmin=0.2):
    ran_cat = fitsio.read(bg_file)
    ran_cat = ran_cat[ran_cat['Z']>zmin]
    # Do the reddening estimate in redshift slices.
    nbins = 5
    est = np.zeros(ran_cat.size)
    est_weight = np.zeros(ran_cat.size)
    zbins = np.percentile(ran_cat['Z'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.

    for i in range(nbins):
        these = (ran_cat['Z'] > zbins[i]) & (ran_cat['Z'] <= zbins[i+1])
        these_est = (Cat.zz > zbins[i]) & (Cat.zz <=zbins[i+1])
        ind = np.random.choice(np.arange(np.sum(these_est)),np.sum(these))
        est[these] = Cat.k[these_est][ind]
        est_weight[these] = Cat.w[these_est][ind]
        
    catalog = treecorr.Catalog(ra=ran_cat['RA'],dec=ran_cat['DEC'],k=est,w=est_weight,\
                                   ra_units='deg',dec_units='deg')
    return catalog
    
def plotres(dd_out=None,fr_out=None, outplotn='fig.png'):
    # Now make a plot.
    dk = fitsio.read(dd_out)
    fr = fitsio.read(fr_out)
    r = np.logspace(-2,2.7,10)
    av = 2.4e-3 * (r/scl)**(-0.84)
    
    fig = plt.figure(figsize=(14,7))
    ### in log space
    ax=fig.add_subplot(121)
    try:
        ax.errorbar(dk['meanr'],dk['kappa'],yerr=dk['sigma'],label='raw')        
        ax.errorbar(dk['meanr'],dk['kappa']-fr['kappa'],yerr=dk['sigma'],label='fr sub')
    except:
        ax.errorbar(dk['meanR'],dk['kappa'],yerr=dk['sigma'],label='raw')        
        ax.errorbar(dk['meanR'],dk['kappa']-fr['kappa'],yerr=dk['sigma'],label='fr sub')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.1,250)
    ax.plot(r,av,label='scaled Menard (2010)') 
    ax.set_xlabel('impact parameter (arcmin)')
    ax.set_ylabel('A_v (mag)')
    ax.legend()    
    ### in linear space
    ax2=fig.add_subplot(122)
    try:
        ax2.errorbar(dk['meanr'],dk['kappa'],yerr=dk['sigma'],label='raw')
        ax2.errorbar(dk['meanr'],dk['kappa'] -fr['kappa'],yerr=dk['sigma'],label='fr sub')
    except:
        ax2.errorbar(dk['meanR'],dk['kappa'],yerr=dk['sigma'],label='raw')        
        ax2.errorbar(dk['meanR'],dk['kappa'] -fr['kappa'],yerr=dk['sigma'],label='fr sub')        
    ax2.plot(r,av,label='scaled Menard (2010)')
    ax2.axhline(0,color='black',linestyle='--',alpha=0.5)
    ax2.set_xlim(0.07,200)
    ax2.set_ylim(-5e-3,0.02)
    ax2.set_xscale('log')
    ax2.legend()
    
    fig.savefig(outplotn)


def get_output_names(basisInd=None,optimal=False):
    if optimal:
        dd_outfile = '../outputs/dustCorr_dd_midz_fg-voptimal.fits'
        dr_outfile = '../outputs/dustCorr_dr_midz_fg-voptimal.fits'
        fr_outfile = '../outputs/dustCorr_fr_midz_fg-voptimal.fits'
        rr_outfile = '../outputs/dustCorr_rr_midz_fg-voptimal.fits'
        fig_outfile = '../outputs/correlationFuncFigures/dustCorr_midz_fg-voptimal.png'       
    elif (basisInd==0):
        dd_outfile = '../outputs/dustCorr_dd_midz_fg-vdust.fits'
        dr_outfile = '../outputs/dustCorr_dr_midz_fg-vdust.fits'
        fr_outfile = '../outputs/dustCorr_fr_midz_fg-vdust.fits'
        rr_outfile = '../outputs/dustCorr_rr_midz_fg-vdust.fits'
        fig_outfile = '../outputs/correlationFuncFigures/dustCorr_midz_fg-vdust.png'      
    else:
        dd_outfile = '../outputs/dustCorr_dd_midz_fg-v'+str(basisInd)+'.fits'
        dr_outfile = '../outputs/dustCorr_dr_midz_fg-v'+str(basisInd)+'.fits'
        fr_outfile = '../outputs/dustCorr_fr_midz_fg-v'+str(basisInd)+'.fits'
        rr_outfile = '../outputs/dustCorr_rr_midz_fg-v'+str(basisInd)+'.fits'
        fig_outfile = '../outputs/correlationFuncFigures/dustCorr_midz_fg-v'+str(basisInd)+'.png'           
    return dd_outfile,dr_outfile,fr_outfile,rr_outfile,fig_outfile



def do_it_all(v,fgCat,fgRan,bgCat,basisInd=None,optimal=False):
    
    dd_outfile,dr_outfile,fr_outfile,rr_outfile,fig_outfile = get_output_names(basisInd,optimal)
    print ("Doing reddening calculation for for vector %s..." % str(v))
    redcat = do_reddening_calculation(bgCat,basisVector=v)
    print ("Done. Getting bg randoms... ")
    print ("Done. Now cross-correlating for vector %s..." % str(v))     
    # Now make the correlation objects.
    DK = treecorr.NKCorrelation(min_sep=0.15,max_sep=200.0,bin_size=.3,sep_units='arcmin')
    DK.process(fgCat,redcat)
    DK.write(dd_outfile)
    FR = treecorr.NKCorrelation(min_sep=0.15,max_sep=200.0,bin_size=.3,sep_units='arcmin')
    FR.process(fgRan,redcat)
    FR.write(fr_outfile)
    plotres(dd_out=dd_outfile,fr_out=fr_outfile,outplotn=fig_outfile)

    return 



def main(argv):
    datapath = '/home/jemcclea/data2/des_dust/catalogs'
    rmz_name = 'DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits'
    rmp_name = 'y1a1-gold-mof-badregion.fits'
    rm_mask = 'DES_Y1A1_3x2pt_redMaGiC_MASK_HPIX4096RING.fits'
    ra_name = 'DES_Y1A1_3x2pt_redMaGiC_RANDOMS.fits'
    #fg_name = 'Sscom_exactArea_galzCut.fits'
    fg_name = 'redMaGiC_midz.fits'
    rmp_file = os.path.join(datapath,rmp_name)
    rmz_file = os.path.join(datapath,rmz_name)
    rmm_file = os.path.join(datapath,rm_mask)
    global ra_file
    ra_file = os.path.join(datapath,ra_name)
    fg_file = os.path.join(datapath,fg_name)
    global zmin
    zmin = 0.5
    global scl
    scl=1.0
    # This parameter decides whether we want to loop over all basis vectors, or use the "optimal" vector
    global optimal
    optimal = False 

    # First, define our orthonormal vector space based on an input extinction vector
    vdust = np.array([1.12224688, 0.82747095, 0.62680647, 0.47880753])
    basis = get_ONbasis(vdust)

    print( "Getting fg catalog and randoms... ")
    fgCat = get_fg_catalog(datapath,fg_file)
    fgRan = get_fg_randoms(maskfile = rmm_file)
    print( "Getting bg science catalog")
    bgCat = get_bg_catalog2(datapath, rmp_file,rmz_file,zmin=zmin)  

    if optimal:      
        print("using optimal dust vector...")
        vec = -0.5923*basis[0]+0.8057175*basis[1]
        
        do_it_all(vec,fgCat,fgRan,bgCat,optimal=True)

    else: 
        # First calculation: "reddening vector"
        vec=basis[0]        
        do_it_all(vec,fgCat,fgRan,bgCat,basisInd=0,optimal=False)
        #Loop through other vectors
        index = 1
        for vec in basis[1:]:
            do_it_all(vec,fgCat,fgRan,bgCat,basisInd=index,optimal=False)
            index+=1

if __name__ == "__main__":
    import pdb, traceback, sys
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)       
