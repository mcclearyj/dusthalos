import numpy as np
from astropy.io import fits
import os
from astropy.table import Table, vstack, hstack
import healpy as hp
import pdb
import fitsio
from esutil import htm
#import astropy.wcs as wcs
import yaml
import re
import os

class AttrDict(dict):
    '''
    More convenient to access dict keys with dict.key than dict['key'],
    so cast the input dict into a class!
    '''

    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def hpRaDecToHEALPixel(ra, dec, nside=4096, nest=False, convert2gal=False):
    phi = ra * np.pi / 180.0
    theta = (90.0 - dec) * np.pi / 180.0

    if convert2gal==True:
        # Transforms ecliptic to galactic coordinates)
        print('Converting ecliptic to galactic')
        r = hp.rotator.Rotator(coord=['E','G'])
        theta_ecl, phi_ecl = r(theta, phi)
        hpInd = hp.ang2pix(nside, theta_ecl, phi_ecl, nest=nest)
    else:
        hpInd = hp.ang2pix(nside, theta, phi, nest=nest)

    return hpInd


def SkyCoordToHEALPixel(skycoord, nside=4096, nest=False):
    '''
    skycoord is an astropy.coordinates.SkyCoord object
    '''
    try:
        phi = skycoord.ra.radian
        theta = np.pi/2. - skycoord.dec.radian
    except:
        phi = skycoord.l.radian
        theta = skycoord.pi/2. - skycoord.b.radian

    hpInd = hp.ang2pix(nside, theta, phi, nest=nest)

    return hpInd



def read_yaml(yaml_file):
    '''
    current package has a problem reading scientific notation as
    floats; see
    https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
    '''
    loader = yaml.SafeLoader
    loader.add_implicit_resolver(
        u'tag:yaml.org,2002:float',
        re.compile(u'''^(?:
        [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$''', re.X),
        list(u'-+0123456789.'))

    with open(yaml_file, 'r') as stream:
        # return yaml.safe_load(stream) # see above issue
        return yaml.load(stream, Loader=loader)


def match_coords(cat1, cat2, ratag1=None, dectag1=None,
                ratag2=None, dectag2=None, radius=0.5):
    '''
    Utility function to match cat1 to cat 2 using celestial coordinates.
    This assumes cat1/cat2 are astropy.Table objects.
    '''

    # Either 'ra/dec' or 'ALPHAWIN_J2000/DELTAWIN_J2000'!
    try:
        if (ratag1 is not None) and (dectag1 is not None):
            cat1_ra = cat1[ratag1]
            cat1_dec =  cat1[dectag1]
        elif 'ra' in cat1.colnames:
            cat1_ra = cat1['ra']
            cat1_dec =  cat1['dec']
        elif 'ALPHAWIN_J2000' in cat1.colnames:
            cat1_ra = cat1['ALPHAWIN_J2000']
            cat1_dec =  cat1['DELTAWIN_J2000']
        else:
            raise KeyError('cat1: no "ra,dec" or ',
                           '"{ALPHA,DELTA}WIN_J2000" columns')
    except:
        raise NameError("\nCouldn't load catalog 1 RA & Dec\n")

    try:
        if (ratag2 is not None) and (dectag2 is not None):
            cat2_ra = cat2[ratag2]
            cat2_dec =  cat2[dectag2]
        elif 'ra' in cat2.colnames:
            cat2_ra = cat2['ra']
            cat2_dec =  cat2['dec']
        elif 'ALPHAWIN_J2000' in cat2.colnames:
            cat2_ra = cat2['ALPHAWIN_J2000']
            cat2_dec =  cat2['DELTAWIN_J2000']
        else:
            raise KeyError('cat2: no "ra,dec" or ',
                           '"{ALPHA,DELTA}WIN_J2000" columns')
    except:
        raise NameError("\nCouldn't load catalog 2 RA & Dec\n")

    cat1_matcher = htm.Matcher(16, ra=cat1_ra, dec=cat1_dec)

    cat2_ind, cat1_ind, dist = cat1_matcher.match(ra=cat2_ra,
                                                  dec=cat2_dec,
                                                  maxmatch=1,
                                                  radius=radius/3600.
                                                  )

    print(f'{len(dist)}/{len(cat1)} gals matched to truth')

    return cat1[cat1_ind], cat2[cat2_ind]
