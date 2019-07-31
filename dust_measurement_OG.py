import numpy as np
import matplotlib.pyplot as plt
import fitsio
import treecorr
import os
import pandas as pd

def get_ortho(vec,index = 0):
    # Get an orthogonal vector with the same norm as that supplied.
    v2 = np.zeros_like(vec)
    v2[index] = 1.0
    v_perp = v2 - np.dot(vec,v2)/np.dot(vec,vec) * vec
    return v_perp*np.sqrt(np.linalg.norm(vec)/np.linalg.norm(v_perp))

def est_reddening(catalog,zeropoint = 24.5, ortho=False,ortho_index = 0):
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
        print "Dust vector is:",vdust
        print "ortho dust vector is:",dmdp
        print "overlap:",np.dot(vdust,dmdp)
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
    return est,wt


def get_fg_catalog(fg_file):
    data = pd.read_csv(fg_file)
    keep = data['bCalCorr'] > 15.0
    catalog = treecorr.Catalog(ra=data['ra'][keep],dec=data['dec'][keep],ra_units='deg',dec_units='deg')

    return catalog

def get_bg_catalog(phot_file,rmz_file,zmin=0.1,ortho=False):
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
    for i in xrange(nbins):
        these = (zcat['ZREDMAGIC'] > zbins[i]) & (zcat['ZREDMAGIC'] <= zbins[i+1])
        this_est,this_wt = est_reddening(cat[these],ortho=ortho)
        est[these] = this_est
        est_weight[these] = this_wt

    catalog = treecorr.Catalog(ra=cat['ra'],dec=cat['dec'],k=est,\
                                   ra_units='deg',dec_units='deg',w=est_weight)
    catalog.zz = zcat['ZREDMAGIC']
    return catalog

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
    for i in xrange(nbins):
        these = (ran_cat['Z'] > zbins[i]) & (ran_cat['Z'] <= zbins[i+1])
        these_est = (Cat.zz > zbins[i]) & (Cat.zz <=zbins[i+1])
        est[these] = np.random.choice(Cat.k[these_est],np.sum(these))
        est_weight[these] = np.random.choice(Cat.w[these_est],np.sum(these))

    catalog = treecorr.Catalog(ra=ran_cat['RA'],dec=ran_cat['DEC'],k=est,w=est_weight,\
                                   ra_units='deg',dec_units='deg')
    return catalog
    


def main(argv):
    datapath = '/home/jemcclea/data2/des_dust'
    rmz_name = 'DES_Y1A1_3x2pt_redMaGiC_zerr_CATALOG.fits'
    rmp_name = 'y1a1-gold-mof-badregion.fits'
    ra_name = 'DES_Y1A1_3x2pt_redMaGiC_RANDOMS.fits'
    fg_name = 'filtered_allsky.csv.gz'
    rmp_file = os.path.join(datapath,rmp_name)
    rmz_file = os.path.join(datapath,rmz_name)
    ra_file = os.path.join(datapath,ra_name)
    fg_file = os.path.join(datapath,fg_name)

    ortho = False
    if ortho:
        index = 0
        dd_outfile = 'dust_correlation_unsubtracted_faint_ortho-'+str(index)+'.fits'
        dr_outfile = 'dust_correlation_randoms_faint_ortho-'+str(index)+'.fits'
    else:
        dd_outfile = 'dust_correlation_unsubtracted_faint.fits'
        dr_outfile = 'dust_correlation_randoms_faint.fits'


    # built a background catalog.
    zmin = 0.15
    print "Getting bg science catalog..."
    bgCat = get_bg_catalog(rmp_file,rmz_file,zmin=zmin,ortho=False)
    print "Done. Getting bg randoms... "
    bgRan = get_bg_randoms(ra_file, bgCat,zmin=zmin)
    print "Done. Getting fg catalog... "
    # Build a foreground catalog.
    fgCat = get_fg_catalog(fg_file)
    print "Done. Now cross-correlating..."
    
    # Now make the correlation objects.
    DK = treecorr.NKCorrelation(min_sep=.5,max_sep=50.0,bin_size=0.2,sep_units='arcmin')
    DK.process(fgCat,bgCat)
    DK.write(dd_outfile)
    RK = treecorr.NKCorrelation(min_sep=.5,max_sep=50.0,bin_size=0.2,sep_units='arcmin')
    RK.process(fgCat,bgRan)
    RK.write(dr_outfile)

    # Now make a plot.
    dk = fitsio.read(dd_outfile)
    dr = fitsio.read(dr_outfile)
    fig,ax = plt.subplots(figsize=(7,7))
    ax.errorbar(dk['meanR'],dk['kappa']-(dr['kappa']-np.mean(dr['kappa'])),yerr=dk['sigma'])
    r = np.logspace(-1,1,10)
    av = 2.5e-3 * (r/2.)**(-0.86)
    ax.plot(r,av,label='Menard (2010)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(1e-4,1)
    ax.set_xlim(.1,100)
    ax.axhline(0,color='black',linestyle='--',alpha=0.5)
    ax.set_xlabel('target separation (arcmin)')
    ax.set_ylabel('A_v (mag)')
    if ortho:
        fig.savefig('dust_extinction-ortho'+str(index)+'.png')
    else:
        fig.savefig('dust_extinction.png')
    
    
if __name__ == "__main__":
    import pdb, traceback, sys
    try:
        main(sys.argv)
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)       
