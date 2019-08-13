import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import fitsio
import treecorr
import os
from astropy.table import Table
import pandas as pd
import healpy as hp


def est_reddening(catalog,zeropoint = 30.0):
    # Make the colors.
    
    gmag = zeropoint - 2.5*np.log10(catalog['mof_flux_g'])
    rmag = zeropoint - 2.5*np.log10(catalog['mof_flux_r'])
    imag = zeropoint - 2.5*np.log10(catalog['mof_flux_i'])
    zmag = zeropoint - 2.5*np.log10(catalog['mof_flux_z'])
    

    data = np.vstack([gmag,rmag,imag,zmag])
    covar = np.cov(data)
    dmdp = np.array([1.12224688, 0.82747095, 0.62680647, 0.47880753])
        
    delta = (np.zeros_like(data).T + dmdp)
    Cinv = np.linalg.inv(covar)
    colors = (data.T - np.average(data,axis=1)).T
    est = np.sum(delta.T*np.dot(Cinv,colors),axis=0)/(np.dot(dmdp,np.dot(Cinv,dmdp)))
    wt = 1./(np.dot(dmdp,np.dot(Cinv,dmdp)))
    
    return est,wt


def get_bg_catalog(phot_file,rmz_file,zmin=0.1):
    joint = os.path.join(datapath,'y1a1_mof_rmz.fits')
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
        cat=fitsio.read(joint,format='fits')
    
    print "background redMaGiC catalog acquired"
    
    # For reference, make color vs. redshift plots of our sample
    print "making redMaGic color vs. redshift plots"
    make_rm_colorVsZ(cat)
    
    # Do the reddening estimate in redshift slices.
    nbins = 10
    
    try:
        est = np.zeros(cat.size)
        est_weight =  np.zeros(cat.size)
    except:
        est = np.zeros(np.shape(cat))
        est_weight =  np.zeros(np.shape(cat))
      
    
    
    zbins = np.percentile(cat['ZREDMAGIC'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.

    for i in range(nbins):
        these = (cat['ZREDMAGIC'] > zbins[i]) & (cat['ZREDMAGIC'] <= zbins[i+1])
        this_est,this_wt = est_reddening(cat[these])
        print ("in slice %d, shape of Av estimator is %d and variance is %f" % (i,len(this_est),this_wt*100.0))
        zrange=[zbins[i],zbins[i+1]]
        make_av_hist(this_est,this_wt*100,zrange)
        make_rm_colormag(cat[these],zrange)
        est[these] = this_est
        est_weight[these] = this_wt

    # Not writing a separate function for this, but plot Av of whole distribution
    try:
        fig=plt.figure(figsize=(8,7))
        ax=fig.add_subplot(111)
        _=ax.hist(est,bins=500,range=[-1.5,1.5])
        ax.set_xlabel('A_v of redMaGic galaxy sample')
        wt_for_plot=np.mean(est_weight) * 100.0
        ax.annotate(('variance ~ %f percent'% wt_for_plot),[.6,.85],xycoords='figure fraction')
        ax.axvline(0,color='black',linestyle='--',alpha=0.5)
        outn='av_redmagicAll.png'
        fig.savefig(outn)
    except:
        pdb.set_trace()
 
    return 0


def make_av_hist(est,wt,zrange):
    fig=plt.figure(figsize=(8,7))
    ax=fig.add_subplot(111)
    _=ax.hist(est,bins=500,range=[-1.5,1.5])
    ax.set_xlabel('A_v of redMagic galaxies at %.3f<z<%.3f'%(zrange[0],zrange[1]))
    ax.annotate(('variance ~ %f percent'% wt),[.6,.85],xycoords='figure fraction')
    ax.axvline(0,color='black',linestyle='--',alpha=0.5)
    outn='av_redmagic_z'+str(np.mean(zrange))+'.png'
    fig.savefig(outn)   
    return 0

def make_rm_colormag(photcat,zrange):
    zeropoint = 30.0
    gmag = zeropoint - 2.5*np.log10(photcat['mof_flux_g'])
    rmag = zeropoint - 2.5*np.log10(photcat['mof_flux_r'])
    imag = zeropoint - 2.5*np.log10(photcat['mof_flux_i'])  
    fig=plt.figure(figsize=(14,7))   
    ax=fig.add_subplot(121)
    ax.plot(rmag,gmag-rmag,',')
    var_gmr=np.std((gmag-rmag))
    ann=ax.annotate("variance is %f.3" % var_gmr)
    ax.set_xlabel('r')
    ax.set_ylabel('g - r')
    ax.set_ylim(1,2.5)
    ax2=fig.add_subplot(122)
    ax2.plot(imag,rmag-imag,',')
    var_rmi=np.std((rmag-imag))
    ann=ax.annotate(("variance is %f.3" % var_rmi))

    ax2.set_xlabel('i')
    ax2.set_ylabel('r - i')
    ax2.set_ylim(0.25,1.75)
    ax2.set_title('redMaGic galaxies at %.3f<z<%.3f'%(zrange[0],zrange[1]))
    outn='redmagic_colormag'+str(np.mean(zrange))+'.png'
    fig.savefig(outn)    
    return 0

def make_rm_colorVsZ(photcat):
    zeropoint = 30.0
    gmag = zeropoint - 2.5*np.log10(photcat['mof_flux_g'])
    rmag = zeropoint - 2.5*np.log10(photcat['mof_flux_r'])
    imag = zeropoint - 2.5*np.log10(photcat['mof_flux_i'])    
    fig=plt.figure(figsize=(14,7))   
    ax=fig.add_subplot(121)
    ax.plot(photcat['ZREDMAGIC'],gmag-rmag,',')
    ax.set_xlabel('redMaGiC redshift')
    ax.set_ylabel('g - r')
    ax.set_ylim(0,3)
    ax2=fig.add_subplot(122)
    ax2.plot(photcat['ZREDMAGIC'],rmag-imag,',')
    ax2.set_xlabel('redMaGiC redshift')
    ax2.set_ylabel('r - i')   
    fig.savefig('redmagic_zVsColor.png')
    
    return 0


def main(argv):
    global datapath
    datapath = '/home/jemcclea/data2/des_dust'
    rmz_name = 'DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits'
    rmp_name = 'y1a1-gold-mof-badregion.fits'
    rm_mask = 'DES_Y1A1_3x2pt_redMaGiC_MASK_HPIX4096RING.fits'
    ra_name = 'DES_Y1A1_3x2pt_redMaGiC_RANDOMS.fits'

    rmp_file = os.path.join(datapath,rmp_name)
    rmz_file = os.path.join(datapath,rmz_name)
    ra_file = os.path.join(datapath,ra_name)

        
    # build a background catalog.
    zmin = 0.15
    print( "Getting bg science catalog...")
    bgCat = get_bg_catalog(rmp_file,rmz_file,zmin=zmin)
    #print ("Done. Getting bg randoms... ")
    #bgRan = get_bg_randoms(ra_file, bgCat,zmin=zmin)
    
    
if __name__ == "__main__":
    import pdb, traceback, sys
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)       
