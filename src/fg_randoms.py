import sys, os
from astropy.table import Table, vstack, hstack
import healpy as hp
import numpy as np
import pdb
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from datetime import datetime
import time

# Local imports
from . import utils
from .hpmask import HpMask

class FgRandoms(HpMask):
        '''
        Subclass of HpMask that produces random galaxies on the sphere. Bit of
        a weird hybrid between HpMask and Catalog.

        Extra attributes
            seed: seed for random number generator; if not set, uses time()
            mask_config: should be able to build a mask from config?
        '''
        def __init__(self, seed=None, config=None, filepath=None,
                    partial=False, coordframe=None):
            self.seed = seed
            self.config = config
            self.rng = None
            self.nrand = 1e6
            self.rand_coords = {}

            # Config check should be run run during runner.
            if config is not None:
                filepath = os.path.join(config['path'], config['filename'])
                coordframe = config['coordframe']

            super().__init__(filepath=filepath, coordframe=coordframe)

        def _generate_rng(self, seed):
            '''
            Quick utility function to generate an rng based on the seed in self,
            though one can be passed too, I guess? Maybe you want to generate
            multiple foregrounds?
            '''

            if (seed == None):
                if self.seed == None:
                    # Didn't set a seed? Set one with time()
                    seed = time.time_ns()
                    self.seed = seed
                    print("\n No seed argument passed, using time.time_ns()\n")
                else:
                    # Use the seed attribute
                    seed = self.seed
                    print(" Using self.seed")
            else:
                print(" Using seed passed to make_fg_randoms()")

            self.rng = np.random.default_rng(seed)
            print(f" Set RNG with seed = {seed}\n")

        def _grab_random_coords(self, nrand):
            '''
            Method to draw random points on the sphere
            '''
            # Did someone set a mask?

            if nrand == None:
                nrand = self.nrand

            # There might actually be a smarter way to do this...
            fcover = np.sum(self.mask > 0)*1./self.mask.size
            ndraw = np.ceil((nrand/fcover)*1.2).astype(int)

            # Generate random pixel indices for the highest resolution
            random_hmap_pixels = self.rng.integers(0,
                                    hp.nside2npix(self.NSIDE), 2*ndraw)

            # Get the RA and Dec coordinates of the random pixels for NSIDE
            rand_ra, rand_dec = hp.pix2ang(self.NSIDE, random_hmap_pixels,
                                    lonlat=True, nest=False)

            self.rand_coords = coord.SkyCoord(ra=rand_ra * u.deg,
                                    dec=rand_dec * u.deg, frame='icrs')

        def write_outfile(self, prefix, mask2=None, overwrite=False):
            '''
            Bit ad-hoc, but write to file
            '''

            # Defile output file name & save
            if self.config is not None:
                outdir = self.config['output_path']
            else:
                outdir = os.path.dirname(self.filepath)
            outcat_name = f'{prefix}{os.path.basename(self.filepath)}'
            outcat_path = os.path.join(outdir, outcat_name)

            # Lazy initialize data table too
            data = Table(self.rand_coords.to_table())
            data.meta['comments'] = []

            # Save specific metadata for single or double-masked data
            if mask2 is not None:
                data.meta['comments'].append(\
                    f'Applied HEALPix masks {self.filepath} and ' + \
                    f'{mask2.filepath} on date {datetime.now():%D %H:%M:%S}')
            else:
                data.meta['comments'].append(
                        f'Applied HEALPix mask {self.filepath} on ' + \
                        f' date {datetime.now():%D %H:%M:%S}')

            data.write(outcat_path, format='fits', overwrite=overwrite)

            print(f'\tSaved catalog to {outcat_path}')


        def make_fg_randoms(self, seed=None, nrand=None):
            '''
            Method to generate random galaxies with masking taken into account.
            '''

            # First, generate the rng
            self._generate_rng(seed=seed)

            # Generate random coodinates with ~appx. the right number of galaxies
            self._grab_random_coords(nrand=nrand)

            # Then write to file
            self.write_outfile(prefix='FgMask_randoms_', overwrite=True)


        def do_overlapping_mask(self, mask2, overwrite=False):
            '''
            Wrapper for apply_overlapping_masks with a mask fed in by a
            runner script. Instantiate background mask, run, save to file.
            '''

            # Get indicies of intersecting coordinates
            good_gals = self.apply_overlapping_masks(coords=self.rand_coords,
                            mask1=self, mask2=mask2)

            # Filter own coordinates
            self.rand_coords = self.rand_coords[good_gals]

            # Write masked to file
            self.write_outfile(prefix='Bg+FgMask_randoms_', mask2=mask2,
                                overwrite=overwrite)
