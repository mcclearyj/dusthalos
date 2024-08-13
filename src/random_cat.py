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
import matplotlib.pyplot as plt

# Local imports
from . import utils
from .hpmask import HpMask

class RandomCat(HpMask):
        """
        Subclass of HpMask that produces random galaxies on the sphere given an
        input HEALPix mask. This is something of a hybrid between HpMask
        and Catalog classes.

        HpMask attributes:
            filepath:  Filepath (relative or absolute) to HEALpix mask
            partial:  Is mask partial? (boolean; default=False)
            NSIDE:  HEALPix NSIDE parameter
            mask:  HpMask map instance
            mask_header:  Mask header info
            seen:  Is HEALPixel filled? (1 or 0; default=0 for False)
            all_nside_hpix:  HEALPix mask, stored as 1 x (12 * NSIDE**2) array
            coordframe: Celestial reference frame (should be recognized
                        astropy.coord.SkyCoord kw like icrs, galactic, ...)

        Extra attributes
            seed: Seed for random number generator; if not set, uses time()
            mask_config: To build a mask from configuration file
            rand_coords: Randomly generated coordinates with mask applied

        """

        def __init__(self, seed=None, mask_config=None, mask_filepath=None,
                        partial=False, mask_coordframe=None, nrand=1e6,
                        output_name=None):
            self.seed = seed
            self.config = mask_config
            self.nrand = nrand
            self.output_name = output_name
            self.rand_coords = {}

            # Config check should be run during runner.
            if mask_config != None:
                mask_filepath = os.path.join(
                    mask_config['path'], mask_config['filename']
                )
                mask_coordframe = mask_config['coordframe']

            super().__init__(filepath=mask_filepath, coordframe=mask_coordframe)

        def _generate_rng(self, seed):
            """
            Quick utility function to generate an rng based on the seed in self,
            though one can be passed too, I guess? Maybe you want to generate
            multiple foregrounds?
            """

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

        def _draw_random_coords(self, nrand):
            """
            Method to draw random points on the sphere by creating catalog of
            uniformly distributed random ras, decs and transforming to mask's
            SkyCoord-format coordframe ('icrs', 'galactic', etc.)

            NOTA BENE: random declinations should be uniform in sin(dec),
                       *not* uniform in dec (unless you like pain).
            """

            # Allow for override of config value
            if nrand == None:
                nrand = self.nrand

            # There might be a smarter way to do this...
            fcover = np.sum(self.mask > 0)*1./self.mask.size
            ndraw = np.ceil((nrand/fcover)*1.2).astype(int)

            """
            # The algorithm below doesn't work!!! If NSIDE is too small, and
            # the HEALPixels are correspondingly large, you end up with large
            # gaps between random points on the sphere, as points will be placed
            # at the center of the HEALPixel area, not randomly within it.

            # Generate random pixel indices for the highest resolution
            random_hmap_pixels = self.rng.integers(
                0, hp.nside2npix(self.NSIDE), 2*ndraw
            )

            # Get the RA and Dec coordinates of the random pixels for NSIDE
            rand_ra, rand_dec = hp.pix2ang(
                self.NSIDE, random_hmap_pixels, lonlat=True, nest=False
            )
            """

            # Sample RA, Dec uniformly on the sphere
            rand_ras = self.rng.uniform(0, 2*np.pi, size=ndraw)
            rand_sindecs = self.rng.uniform(
                np.sin(-np.pi/2), np.sin(np.pi/2), size=ndraw
            )
            rand_decs = np.arcsin(rand_sindecs)
            rand_ras = np.rad2deg(rand_ras); rand_decs = np.rad2deg(rand_decs)

            # Create SkyCoord instance from random RA/Dec pairs above
            icrs_rand_coords = coord.SkyCoord(
                ra=rand_ras * u.deg, dec=rand_decs * u.deg, frame='icrs'
            )

            # Populate rand_coords in specified celestial coordinate frame
            self.rand_coords = icrs_rand_coords.transform_to(self.coordframe)

        def plot_radec(self, plot_config):
            """
            FOR DEBUGGING PURPOSES: plot RA/Dec distributions.
            """

            # We are comparing to another, outside catalog
            if catalog_for_comparison != None:
                comparison_cat = Table.read(catalog_for_comparison, memmap=True)
                rand_ras = comparison_cat[comparison_ra_key]
                rand_decs = comparison_cat[comparison_dec_key]
                cc_label = f'{os.path.basename(catalog_for_comparison)}'

            # Just make uniform points on sphere
            else:
                rand_ras = self.rng.uniform(0, 2*np.pi, size=int(self.nrand))
                rand_ras = np.rad2deg(rand_ras)
                rand_sindecs = self.rng.uniform(
                    np.sin(-np.pi/2), np.sin(np.pi/2), size=int(self.nrand)
                )
                rand_decs = np.rad2deg(np.arcsin(rand_sindecs))
                cc_label = 'Uniform distribution'

            # Create coord bins for histogram plots (unit: degrees)
            nn, ra_bins = np.histogram(self.rand_coords.icrs.ra.deg, bins=200)
            nn, dec_bins = np.histogram(self.rand_coords.icrs.dec.deg, bins=200)

            # Create figure instance
            fig, axs = plt.subplots(1, 2, figsize=(14,6), tight_layout=True)

            # Plot right ascensions
            axs[0].hist(
                rand_ras, bins=ra_bins, alpha=0.7, histtype='stepfilled',
                density=True, label=cc_label
            )
            axs[0].hist(
                self.rand_coords.icrs.ra.deg, bins=ra_bins, alpha=0.8,
                histtype='stepfilled', density=True, label='masked random cat'
            )
            axs[0].set_xlabel('RA'); axs[0].set_ylabel('Probability Density')
            axs[0].legend(loc='lower right')

            # Plot declinations
            axs[1].hist(
                rand_decs, bins=dec_bins, alpha=0.7, histtype='stepfilled',
                density=True, label=cc_label
            )
            axs[1].hist(
                self.rand_coords.icrs.dec.deg, bins=dec_bins, alpha=0.8,
                histtype='stepfilled', density=True, label='random cat'
            )
            axs[1].set_xlabel('Dec'); axs[0].set_ylabel('Probability Density')
            axs[1].legend(loc='lower right')

            # Save plot to file
            fig.savefig(os.path.join(self.config['path'], plotname))

        def write_outfile(self, prefix, mask2=None, overwrite=False):
            """
            Write random coordinates to file
            """

            # Define output file name & directory
            if self.output_name != None:
                outcat_name = f'{prefix}{self.output_name}'
            else:
                outcat_name = f'{prefix}{os.path.basename(self.filepath)}'

            outdir = os.path.dirname(self.filepath)
            outcat_path = os.path.join(outdir, outcat_name)

            # Lazy-initialize data table to store coordinates
            data = Table(self.rand_coords.to_table())
            data.meta['comments'] = []

            # Save metadata for single or double-masked data, write to file
            if mask2 != None:
                data.meta['comments'].append(\
                    f'Applied HEALPix masks {self.filepath} and ' + \
                    f'{mask2.filepath} on date {datetime.now():%D %H:%M:%S}')
            else:
                data.meta['comments'].append(
                        f'Applied HEALPix mask {self.filepath} on ' + \
                        f' date {datetime.now():%D %H:%M:%S}')

            data.write(outcat_path, format='fits', overwrite=overwrite)

            print(f'\tSaved catalog to {outcat_path}')

        def make_random_cat(self, seed=None, nrand=None):
            """
            Method to generate random galaxies with masking taken into account.
            """

            # First, generate the rng
            self._generate_rng(seed=seed)

            # Generate random coodinates with ~appx. the right number of galaxies
            self._draw_random_coords(nrand=nrand)

            # Apply own mask!
            seen = self.apply_mask(coords=self.rand_coords)
            self.rand_coords = self.rand_coords[seen]

            # Then write to file
            self.write_outfile(prefix='', overwrite=True)

        def do_overlapping_mask(self, mask2, overwrite=False,
                                catalog_for_comparison=None):
            """
            Wrapper for apply_overlapping_masks with a mask fed in by a
            runner script. Instantiate background mask, run, save to file.
            """

            # Get indicies of intersecting coordinates
            good_gals = self.apply_overlapping_masks(
                coords=self.rand_coords, mask1=self, mask2=mask2
            )

            # Filter own coordinates
            self.rand_coords = self.rand_coords[good_gals]

            # Plotname isn't an argument right now; but maybe in the future
            if catalog_for_comparison != None:
                self.plot_radec(
                    plotname='masked_randoms_RADEC_distrib_comparison.png',
                    catalog_for_comparison=catalog_for_comparison,
                    comparison_ra_key='ra', comparison_dec_key='dec'
                )

            # Write masked to file
            self.write_outfile(
                prefix='BgFgMasked_', mask2=mask2, overwrite=overwrite
            )
