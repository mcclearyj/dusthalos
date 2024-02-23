# Plot fiducial vs. "best-fit" dust vectors against each other

def get_ONbasis(vdust):
    # construct an orthonormal set of basis vectors according
    # to the Gram-Schmidt procedure: define V={v0,v1,v2,v3}
    # that spans the "extinction" basis

    # right now, we have a single vector along which (we think) reddening occurs
    # so define ON basis vectors along which we can project the (real?) vector

    # start by normalizing the input dust vector
    vec = vdust#/norm(vdust)
    
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

def make_figure():
    vdust = np.array([1.12224688, 0.82747095, 0.62680647, 0.47880753])
    basis=get_ONbasis(vdust)

    ## Coefficients obtained by eval_chi2.py on the kfid/k1/k2/k3 obtained by running dust_measurement.py
    ## with the "standard" ON_basis vectors
    new=basis[0]-0.302*basis[1]-0.052*basis[2]-0.14*basis[3]
    new_norm=new/1.23


    # define array of central wavelengths for plotting 
    wl=np.array([475.4,620.4,769.8,966.5])

    fig=plt.figure(figsize=(7,7))
    ax=fig.add_subplot(111)
    ax.plot(wl, fid,'.b',markersize=15,label='fiducial dust vector')
    ax.plot(wl, fid,'--k')
    ax.plot(wl, new_norm, '.r',markersize=15,label='true (?) dust vector')
    ax.plot(wl, new_norm,'--k')

    ax.axvline(537.7)
    ax.axhline(1)

    ax.set_xlabel('wavelength',fontsize=16)
    ax.set_ylabel('A_lambda/A_V',fontsize=16)
    for tick in ax.xaxis.get_major_ticks():  
     tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
     tick.label.set_fontsize(16)

    ax.legend()

    return

