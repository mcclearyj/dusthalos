import sys, os
import healpy as hp
import numpy as np
import pdb
from datetime import datetime
import time
import matplotlib.pyplot as plt

### Astropy imports
from astropy.io import fits
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack, hstack

### For multiprocessing 
import multiprocessing as mp

### Local imports
from . import utils
from .hpmask import HpMask
# ----------------------------------------------------------------
# ---------------- GLOBALS FOR MULTIPROCESSING -------------------
_WORKER = {}  # process-local in each child

def _init_worker(mask_bool, nside, nest, seed_base):
    """Initialize worker process with necessary parameters."""
    global _WORKER
    _WORKER['mask']  = mask_bool
    _WORKER['nside'] = int(nside)
    _WORKER['nest']  = bool(nest)
    seed = (int(seed_base) ^ (os.getpid() & 0xFFFFFFFF)) & 0xFFFFFFFF
    _WORKER['rng']   = np.random.default_rng(seed)

def _draw_and_filter(batch_size):
    rng   = _WORKER['rng']           # <-- same dict the initializer filled
    mask  = _WORKER['mask']
    nside = _WORKER['nside']
    nest  = _WORKER['nest']

    # draw uniformly on sphere (degrees)
    lon = rng.uniform(0.0, 360.0, size=batch_size)
    s   = rng.uniform(-1.0, 1.0, size=batch_size)   # sin(lat)
    lat = np.degrees(np.arcsin(s))

    pix  = hp.ang2pix(nside, lon, lat, lonlat=True, nest=nest)
    keep = mask[pix]
    return lon[keep], lat[keep]

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
            self._worker = {}

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

        # ---------------- uniform-on-sphere draws in degrees ----------------
        def _draw_uniform_lonlat_deg(n, rng):
            """
            Uniform on the sphere.
            Returns lon, lat in DEGREES:
            lon ∈ [0, 360), lat ∈ [-90, 90]
            """
            lon = rng.uniform(0.0, 360.0, size=n)
            s = rng.uniform(-1.0, 1.0, size=n)    # sin(lat)
            lat = np.degrees(np.arcsin(s))
            return lon, lat

        def sample_in_mask_parallel(
                self, 
                nrand, 
                *, 
                seed=None, 
                n_workers=None, 
                batch_per_worker=500_000
        ):
            mask_bool = (self.mask > 0) if self.mask.dtype != bool else self.mask
            nside     = hp.get_nside(mask_bool)
            nest      = getattr(self, "mask_is_nest", False)
            if seed is None:
                seed = int(np.random.SeedSequence().generate_state(1)[0])

            if n_workers is None:
                n_workers = max(1, mp.cpu_count() - 1)

            acc_lon, acc_lat = [], []
            need = int(nrand)

            ctx = mp.get_context("spawn")
            with ctx.Pool(processes=n_workers,
                        initializer=_init_worker,
                        initargs=(mask_bool, nside, nest, seed)) as pool:
                while need > 0:
                    for lon_chunk, lat_chunk in pool.map(_draw_and_filter, [batch_per_worker]*n_workers):
                        take = min(need, lon_chunk.size)
                        if take:
                            acc_lon.append(lon_chunk[:take])
                            acc_lat.append(lat_chunk[:take])
                            need -= take
                        if need == 0:
                            break

            lon_deg = np.concatenate(acc_lon)
            lat_deg = np.concatenate(acc_lat)
            # build SkyCoord once (in the mask’s frame)
            if str(self.coordframe).lower() == 'galactic':
                sc = SkyCoord(l=lon_deg*u.deg, b=lat_deg*u.deg, frame='galactic')
            else:
                sc = SkyCoord(ra=lon_deg*u.deg, dec=lat_deg*u.deg, frame='icrs')
            return sc

        # ---------------- Class-friendly wrapper ----------------
        def _draw_random_coords_parallel(self, nrand=None, seed=None, n_workers=None, batch_per_worker=500_000):
            """
            Draw directly in the mask's frame (self.coordframe), filter via HEALPix,
            and store SkyCoord in that same frame.
            """
            if nrand is None:
                nrand = int(self.nrand)

            # Determine the mask frame string the same way your apply_mask does
            # (you've been using self.coordframe == 'galactic' to branch)
            mask_frame = 'galactic' if str(self.coordframe).lower() == 'galactic' else 'icrs'

            # NESTED or RING? Your apply_mask used nest=False; keep the same unless you know otherwise.
            nest = getattr(self, "mask_is_nest", False)

            self.rand_coords = self.sample_in_mask_parallel(
                nrand=nrand,
                n_workers=n_workers,
                batch_per_worker=batch_per_worker,
                seed=seed,
            )

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

        def plot_radec(self, plotname, catalog_for_comparison=None, comparison_ra_key=None, comparison_dec_key=None):
            """
            FOR DEBUGGING PURPOSES: plot RA/Dec distributions.
            """

            # We are comparing to another, outside catalog
            if catalog_for_comparison is not None:
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
            #self._generate_rng(seed=seed)

            # Generate random coodinates with ~appx. the right number of galaxies
            #self._draw_random_coords(nrand=nrand)
            self._draw_random_coords_parallel(
                nrand=nrand, seed=seed, 
                n_workers=None, batch_per_worker=500_000
            )

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
