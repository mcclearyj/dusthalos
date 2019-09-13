import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pdb
import fitsio
import os

def chi2(alpha,kappa_norm,sigma,kappa1_norm,sigma1,sigma_sys = 1e-4):
      numerator=(alpha*kappa_norm - kappa1_norm)**2
      denominator=(alpha**2*sigma**2 +sigma1**2 + sigma_sys**2)
      x2 = np.sum(np.divide(numerator,denominator))
      return x2

def chi2_beta(beta,kappa_norm,sigma,kappa1_norm,sigma1,sigma_sys = 1e-4):
      """ Now constant of proportionality is on kappa1 """
      numerator=(kappa_norm - beta*kappa1_norm)**2
      denominator=(sigma**2 +beta**2*sigma1**2+sigma_sys**2)
      x2 = np.sum(np.divide(numerator,denominator))
      return x2

def make_a_chi2_plot(parr,chiarr):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(alpha_arr,chi2_arr)
    fig.savefig('chi2_arr.png')   
    return 0

def main():

    # uncomment your galaxy sample of choice
    #sample = "iifsc"
    #sample = "galex"
    sample = "scos"
    print("using %s sample\n" % sample)
    outpath = '../outputs'
    if (sample=='iifsc'):
        
        ddust=fitsio.read(os.path.join(outpath,'dustCorr_dd_orthonorm-vdust-iifsc.fits')) # The dust signal -- real fg x real background 
        frdust=fitsio.read(os.path.join(outpath,'dustCorr_fr_orthonorm-vdust-iifsc.fits')) # random fg x background 
        dd1=fitsio.read(os.path.join(outpath,'dustCorr_dd_orthonorm-v1-iifsc.fits')) # The dust signal -- real fg x real background 
        fr1=fitsio.read(os.path.join(outpath,'dustCorr_fr_orthonorm-v1-iifsc.fits')) # random fg x background
        
    elif(sample=='galex'):
        
        ddust=fitsio.read(os.path.join(outpath,'dustCorr_dd_orthonorm-vdust-galex.fits')) # The dust signal -- real fg x real background 
        frdust=fitsio.read(os.path.join(outpath,'dustCorr_fr_orthonorm-vdust-galex.fits')) # random fg x background
        dd1=fitsio.read(os.path.join(outpath,'dustCorr_dd_orthonorm-v1-galex.fits')) # The dust signal -- real fg x real background 
        fr1=fitsio.read(os.path.join(outpath,'dustCorr_fr_orthonorm-v1-galex.fits')) # random fg x background

    elif(sample=='scos'):
        ## Assuming you want Scos?
        ddust=fitsio.read(os.path.join(outpath,'dust_correlation_dd_orthonorm-vdust.fits')) # The dust signal -- real fg x real background 
        frdust=fitsio.read(os.path.join(outpath,'dust_correlation_fr_orthonorm-vdust.fits')) # random fg x background
        dd1=fitsio.read(os.path.join(outpath,'dust_correlation_dd_orthonorm-v1.fits')) # The dust signal -- real fg x real background 
        fr1=fitsio.read(os.path.join(outpath,'dust_correlation_fr_orthonorm-v1.fits')) # fg x random background

    else:
        print( "Sample not specified?")
        pdb.set_trace()

    # Now make "normalized" i.e. foreground-subtracted kappas
    kappaDust_norm = ddust['kappa']-frdust['kappa']
    sigmaDust=ddust['sigma']
    kappa1_norm = dd1['kappa']-fr1['kappa']
    sigma1=dd1['sigma']

    # generate an alpha parameter array
    print("#\n#minimizing the chi2 function (alpha*kappa_norm - kappa1_norm)**2\n#")
    alpha_arr=np.arange(-2, 1, 0.001)
    chi2_arr=np.zeros_like(alpha_arr)
    chi_arr = []
    for i,alpha in enumerate(alpha_arr):
        chi2_arr[i]=chi2(alpha,kappaDust_norm,sigmaDust,kappa1_norm,sigma1)
    alphabest = alpha_arr[chi2_arr==min(chi2_arr)]
    print("%f is best-fit alpha of propotionality" % alphabest)   
    print("reduced chi2 for alpha is %f\n" %(min(chi2_arr)/12.))
    

    # do same thing, but with a beta on the k1 values
    print("#\n#minimizing the chi2 (beta*kappa_norm - kappa1_norm)**2\n#")
    beta_arr=np.arange(0, 4, 0.001)
    chi2_arr=np.zeros_like(beta_arr)
    for i,beta in enumerate(beta_arr):
        chi2_arr[i]=chi2_beta(beta,kappaDust_norm,sigmaDust,kappa1_norm,sigma1)       
    print("%f is best-fit beta of propotionality" % (beta_arr[chi2_arr==min(chi2_arr)]))   
    print("reduced chi2 for beta is %f\n" %(min(chi2_arr)/12.))
    pdb.set_trace()

    
    
if __name__ == "__main__":
    import pdb, traceback, sys
    try:
        main()
    except:
        thingtype, value, tb = sys.exc_info()
        traceback.print_exc()
        pdb.post_mortem(tb)       
