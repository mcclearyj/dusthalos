import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import fitsio
import treecorr
import os
from astropy.table import Table
import pandas as pd
import healpy as hp


def get_ortho(vec,index = 0):
    # Get an orthogonal vector with the same norm as that supplied.
    v2 = np.zeros_like(vec)
    v2[index] = 1.0
    v_perp = v2 - np.dot(vec,v2)/np.dot(vec,vec) * vec
    return v_perp*np.sqrt(np.linalg.norm(vec)/np.linalg.norm(v_perp))

def est_reddening(catalog,zeropoint = 30.0, ortho=False,ortho_index = 0):
    # Make the colors.
    
    gmag = zeropoint - 2.5*np.log10(catalog['mof_flux_g'])
    rmag = zeropoint - 2.5*np.log10(catalog['mof_flux_r'])
    imag = zeropoint - 2.5*np.log10(catalog['mof_flux_i'])
    zmag = zeropoint - 2.5*np.log10(catalog['mof_flux_z'])
    

    data = np.vstack([gmag,rmag,imag,zmag])
    covar = np.cov(data)
    if ortho:
        #vdust = np.array([ 4.60, 3.10, 2.18, 1.6 ])
        # Fitzpatrick99
        # vdust = np.array([1.12150099, 0.77164321, 0.57725486, 0.45259124])
        # Calzetti:
        #vdust = np.array([1.13552323, 0.77956032, 0.55082914, 0.39834168])
        # Cardelli, Clayton & Mathis '89
        vdust = np.array([1.12224688, 0.82747095, 0.62680647, 0.47880753])
        dmdp = get_ortho(vdust, index = ortho_index)
        print( "Dust vector is:",vdust)
        print ("ortho dust vector is:",dmdp)
        print ("overlap:",np.dot(vdust,dmdp))
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
    est = np.sum(delta.T*np.dot(Cinv,colors),axis=0)/(np.dot(dmdp,np.dot(Cinv,dmdp)))
    wt = 1./(np.dot(dmdp,np.dot(Cinv,dmdp)))
    
    return est,wt


def get_bg_catalog(phot_file,rmz_file,zmin=0.1,ortho=False,ortho_index=-1):
    bg_phot = fitsio.read(phot_file)
    bg_rmz = fitsio.read(rmz_file)

    # Match these.
    inds = np.in1d(bg_phot['coadd_objects_id'],bg_rmz['ID'])
    cat = bg_phot[inds]
    cat.sort(order='coadd_objects_id')
    bg_rmz.sort(order='ID')
    
    
    # Filter these.
    keep = (cat['mof_flags']==0) & (cat['flags_badregion']==0) & \
           (cat['mof_flux_g']>0) & (cat['mof_flux_r']>0) & \
           (cat['mof_flux_i']>0) & (cat['mof_flux_z']>0) & \
           (bg_rmz['ZREDMAGIC']>zmin)
    cat = cat[keep]
    zcat = bg_rmz[keep]
       
    # Do the reddening estimate in redshift slices.
    nbins = 10
    est = np.zeros(cat.size)
    est_weight = np.zeros(cat.size)
    
    zbins = np.percentile(zcat['ZREDMAGIC'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.

    for i in range(nbins):
        these = (zcat['ZREDMAGIC'] > zbins[i]) & (zcat['ZREDMAGIC'] <= zbins[i+1])
        this_est,this_wt = est_reddening(cat[these],ortho=ortho,ortho_index = ortho_index)
        est[these] = this_est
        est_weight[these] = this_wt

    catalog = treecorr.Catalog(ra=cat['ra'],dec=cat['dec'],k=est,\
                                  ra_units='deg',dec_units='deg',w=est_weight)
    
    catalog.zz = zcat['ZREDMAGIC']

    return catalog

def get_bg_catalog2(datapath,phot_file,rmz_file,zmin=0.15,ortho=False,ortho_index = -1):

    """ Use when the matching catalog has already been made.
        TO DO: merge with function above, make it a try/except
    """
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

    print "generating bg catalog for correlation"

    # Do the reddening estimate in redshift slices.                                                                                                                                          
    nbins = 10
    est = np.zeros(cat.size)
    est_weight = np.zeros(cat.size)

    zbins = np.percentile(cat['ZREDMAGIC'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.

    for i in range(nbins):
        these = (cat['ZREDMAGIC'] > zbins[i]) & (cat['ZREDMAGIC'] <= zbins[i+1])
        this_est,this_wt = est_reddening(cat[these],ortho=ortho,ortho_index = ortho_index)
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


def get_fg_catalog(fg_filen,maskfile = None,nside=4096,nest=False):
    
    #data = Table.read(fg_filen,format='fits')
    #ra_cat=data['RAJ2000']
    #dec_cat=data['DEJ2000']
    
    try:
        data = Table.read(fg_filen,format='csv')

        # To trim catalog down significantly, 
        # and select on redshift, cut on magnitude
        wg,=np.where((data['nuv_mag'] <= 21.5) & (data['nuv_class_star'] <= 0.5) )
        data=data[wg]

        ra_cat =  data['ra']
        dec_cat = data['dec']

        # First, determine indices of maskfile that are non-zero
        hmap = hp.read_map(maskfile,nest=False)
        keep = hmap != hp.UNSEEN
        hpIndex=np.arange(np.size(hmap))
    
        catHpInd = hpRaDecToHEALPixel(ra_cat,dec_cat,nside=nside,nest=nest)
        
        fgKeep=np.in1d(catHpInd,hpIndex[keep])

        # If desired, write this trimmed GALEX catalog to file
        data[fgKeep].write('galex_trimmed.fits',format='fits',overwrite=True)
        tcCatalog = treecorr.Catalog(ra=data[fgKeep]['ra'],dec=data[fgKeep]['dec'],\
                                                   ra_units='deg',dec_units='deg')
        print( "Length of catalog after cuts = %i" % len(fgKeep.nonzero()[0]))
    except:
        # Hopefully this is the trimmed catalog!
        data = Table.read(fg_filen,format='fits')
        try:
            ra_cat =  data['ra']
            dec_cat = data['dec']
        except:
            ra_cat = data['_RAJ2000']
            dec_cat = data['_DEJ2000']
        tcCatalog = treecorr.Catalog(ra=ra_cat,dec=dec_cat,ra_units='deg',dec_units='deg')
        print( "Length of catalog read in = %i" % len(dec_cat))

    return tcCatalog

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
    nbins = 10
    est = np.zeros(ran_cat.size)
    est_weight = np.zeros(ran_cat.size)
    zbins = np.percentile(ran_cat['Z'],np.linspace(0,100,nbins+1))
    zbins[0] = 0.
    zbins[-1] = zbins[-1] + 1.
    #zbins=[i*0.15 for i in range(6)];zbins=[zbins[i] + .15 for i in range(len(zbins))]

    for i in range(nbins):
        these = (ran_cat['Z'] > zbins[i]) & (ran_cat['Z'] <= zbins[i+1])
        these_est = (Cat.zz > zbins[i]) & (Cat.zz <=zbins[i+1])
        ind = np.random.choice(np.arange(np.sum(these_est)),np.sum(these))
        est[these] = Cat.k[these_est][ind]
        est_weight[these] = Cat.w[these_est][ind]
        
    catalog = treecorr.Catalog(ra=ran_cat['RA'],dec=ran_cat['DEC'],k=est,w=est_weight,\
                                   ra_units='deg',dec_units='deg')
    return catalog
    
def plotres(dd_out,dr_out,fr_out=None,rr_out = None, ortho = False,ortho_index = 0,scl=1.0):
    # Now make a plot.
    dk = fitsio.read(dd_out)
    dr = fitsio.read(dr_out)
    fr = fitsio.read(fr_out)
    rr = fitsio.read(rr_out)
    
    fig = plt.figure(figsize=(14,7))
    ### in log space
    ax=fig.add_subplot(121)
    try:
        ax.errorbar(dk['meanr'],dk['kappa'],yerr=dk['sigma'],label='raw')        
        #ax.errorbar(dk['meanr'],dk['kappa']-fr['kappa'],yerr=dk['sigma'],label='fg random subtraction')
        #ax.errorbar(dk['meanr'],dk['kappa'] -fr['kappa'] - dr['kappa'] + rr['kappa'],yerr=dk['sigma'],label='LZ++')
    except:
        ax.errorbar(dk['meanR'],dk['kappa']-fr['kappa'],yerr=dk['sigma'],label='raw')        
        #ax.errorbar(dk['meanR'],dk['kappa']-fr['kappa'],yerr=dk['sigma'],label='fg random subtraction')
        #ax.errorbar(dk['meanR'],dk['kappa'] -fr['kappa'] - dr['kappa'] + rr['kappa'],yerr=dk['sigma'],label='LZ++')
    ax.set_xscale('log')
    ax.set_yscale('log')
    #ax.set_ylim(1e-6,.2)
    ax.set_xlim(0.1,250)
    r = np.logspace(-2,2.7,10)
    av = 2.4e-3 * (r/scl)**(-0.84)
    ax.plot(r,av,label='Menard (2010)') 
    ax.set_xlabel('impact parameter (arcmin)')
    ax.set_ylabel('A_v (mag)')
    ax.legend()

    ### in linear space
    ax2=fig.add_subplot(122)
    try:
        ax2.errorbar(dk['meanr'],dk['kappa'],yerr=dk['sigma'],label='raw')
        #ax2.errorbar(dk['meanr'],dk['kappa']-fr['kappa'],yerr=dk['sigma'],label='fg random subtraction')
        ax2.errorbar(dk['meanr'],dk['kappa'] -fr['kappa'] - dr['kappa'] + rr['kappa'],yerr=dk['sigma'],label='LZ++')
    except:
        ax2.errorbar(dk['meanR'],dk['kappa'],yerr=dk['sigma'],label='raw')        
        #ax2.errorbar(dk['meanR'],dk['kappa']- fr['kappa'],yerr=dk['sigma'],label='fg random subtraction')
        ax2.errorbar(dk['meanR'],dk['kappa'] -fr['kappa'] - dr['kappa'] + rr['kappa'],yerr=dk['sigma'],label='LZ++')
        
    ax2.plot(r,av,label='scaled Menard (2010)')
    ax2.axhline(0,color='black',linestyle='--',alpha=0.5)
    ax2.set_xlim(0.07,200)
    ax2.set_ylim(-5e-4,0.06)
    ax2.set_xscale('log')
    ax2.legend()
    if not ortho:
        fig.savefig('../outputs/correlationFuncFigures/dust_corr.png')
    else:
        fig.savefig('../outputs/correlationFuncFigures/iifscDust_corr_ortho-'+str(ortho_index)+'.png')

def main(argv):
    datapath = '/home/jemcclea/data2/des_dust/catalogs'
    rmz_name = 'DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits'
    rmp_name = 'y1a1-gold-mof-badregion.fits'
    rm_mask = 'DES_Y1A1_3x2pt_redMaGiC_MASK_HPIX4096RING.fits'
    ra_name = 'DES_Y1A1_3x2pt_redMaGiC_RANDOMS.fits'
    #fg_name='galex_trimmed.fits'
    fg_name='iifsc_des_overlap.fits'
    scl = 13.2 # This changes depending on avg. redshift of fg
    #scl = 1.92
    #fg_name = 'Sscom_exactArea_galzCut.fits'
    rmp_file = os.path.join(datapath,rmp_name)
    rmz_file = os.path.join(datapath,rmz_name)
    rmm_file = os.path.join(datapath,rm_mask)
    ra_file = os.path.join(datapath,ra_name)
    fg_file = os.path.join(datapath,fg_name)
    plot = True
    global ortho
    ortho = True
    index = 3

    
    if ortho:

        dd_outfile = '../outputs/dust_correlation_dd_ortho-'+str(index)+'-iifsc.fits'
        dr_outfile = '../outputs/dust_correlation_dr_ortho-'+str(index)+'-iifsc.fits'
        fr_outfile = '../outputs/dust_correlation_fr_ortho-'+str(index)+'-iifsc.fits'
        rr_outfile = '../outputs/dust_correlation_rr_ortho-'+str(index)+'-iifsc.fits'
    else:
        dd_outfile = '../outputs/dust_correlation_dd10.fits'
        dr_outfile = '../outputs/dust_correlation_dr10.fits'
        fr_outfile = '../outputs/dust_correlation_fr10.fits'
        rr_outfile = '../outputs/dust_correlation_rr10.fits'        
        
    # build a background catalog.
    zmin = 0.2
    print("Getting fg catalog...")
    fgCat = get_fg_catalog(fg_file,maskfile = rmm_file)

    print( "Done. Getting bg science catalog...")
    #bgCat = get_bg_catalog(rmp_file,rmz_file,zmin=zmin,ortho=ortho,ortho_index = index)
    # use if we have already done redMaGiC/redshift xcorr
    bgCat = get_bg_catalog2(datapath,rmp_file,rmz_file,zmin=zmin, ortho=ortho,ortho_index = index)
    print ("Done. Getting bg randoms... ")
    bgRan = get_bg_randoms(ra_file, bgCat,zmin=zmin)
    print( "Done. Getting fg randoms... ")
    # Build a foreground random catalog, keep only what overlaps the DES coverage
    fgRan = get_fg_randoms(maskfile = rmm_file)
    
    print ("Done. Now cross-correlating...")
    
    # Now make the correlation objects.

    DK = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=.6,sep_units='arcmin')
    DK.process(fgCat,bgCat)
    DK.write(dd_outfile)
    RK = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=.6,sep_units='arcmin')
    RK.process(fgCat,bgRan)
    RK.write(dr_outfile)
    FR = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=.6,sep_units='arcmin')
    FR.process(fgRan,bgCat)
    FR.write(fr_outfile)
    RR = treecorr.NKCorrelation(min_sep=0.1,max_sep=200.0,bin_size=0.6,sep_units='arcmin')
    RR.process(fgRan,bgRan)
    RR.write(rr_outfile)
    

    if plot:
        plotres(dd_outfile,dr_outfile,fr_out = fr_outfile,rr_out=rr_outfile,ortho = ortho,ortho_index=index)

    
if __name__ == "__main__":
    import pdb, traceback, sys
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)       
