import statsmodels.api as sm

def do_WLS_fit(this_dk, dust_rel, name=None, start=None, end=None, log=False):
    """
    Perform a weighted least-squares linear regression to
    data, taking care to transform the variables appropriately.

    Inputs
        this_dk: should be the 'compensated' treecorr output
        dust_rel: the corrected A_V relationship (dk-dr+rr)
        name: what are we calling it?
        start: exclude points with index < 'start' for fitting
        end: ibidem, but for points with index > 'end'
        log: if True, fit the relationship in log-space to stabilize the fits.
             Default is False.
    """
    dust_rel = dust_rel[start:end]
    this_dk = this_dk[start:end]
    ## Transform the dependent variable (fully compensated kappa)
    if log == True:
        # Keep positive values only...
        ok = (dust_rel > 0)
        dust_rel = dust_rel[ok]
        this_dk = this_dk[ok]
        # Fit log of dust relation
        Y = np.log(dust_rel)
        ## Also need to transform the weights!
        weights_OG = 1/(this_dk['sigma']**2)
        weight_scale = (1/Y)**2
        weights = weights_OG * weight_scale

    ## Just do regular fit
    else:
        Y = dust_rel
        weights = 1/(this_dk['sigma']**2)

    ## If plotting X on a log-scale, which we are,
    ## fit relationship on log scale
    X = this_dk['meanlogr']

    ## Add intercept to abscissa
    X = sm.add_constant(X)

    ## Do fit
    mod_wls = sm.WLS(Y, X, weights=weights)
    res_wls = mod_wls.fit()

    ## Print results
    print("")
    print(f"Fit result: {name}")
    print("")
    print(res_wls.summary())

    return(res_wls)
