import os

import numpy as np
from astropy.table import Table

from .hpmask import HpMask


class DemeanTreecorrCat:
    """Load a TreeCorr catalog, subtract the associated HEALPix mean value, and save it."""

    def __init__(
        self,
        cat_dir,
        cat_name,
        map_file,
        nside=256,
        output_name=None,
        ra_col="ra",
        dec_col="dec",
        value_col="k",
        output_dir=None,
        coord_frame="icrs",
    ):
        self.cat_dir = cat_dir
        self.cat_name = cat_name
        self.map_file = map_file
        self.nside = nside
        self.output_name = output_name or f"demeaned_{os.path.splitext(cat_name)[0]}.fits"
        self.ra_col = ra_col
        self.dec_col = dec_col
        self.value_col = value_col
        self.output_dir = output_dir or cat_dir
        self.coord_frame = coord_frame

    @property
    def cat_path(self):
        return os.path.join(self.cat_dir, self.cat_name)

    @property
    def map_path(self):
        return os.path.join(self.cat_dir, self.map_file)

    @property
    def output_path(self):
        return os.path.join(self.output_dir, self.output_name)

    def load_catalog(self):
        return Table.read(self.cat_path, memmap=True, format="fits")

    def demean_catalog(self, catalog=None):
        """Return a copy of the catalog with a de-meaned value column added."""
        if catalog is None:
            catalog = self.load_catalog()

        catalog = catalog.copy()

        ra = catalog[self.ra_col].data
        dec = catalog[self.dec_col].data
        values = catalog[self.value_col].data

        cat_pix = HpMask.coords_to_healpixels(
            lon=ra,
            lat=dec,
            frame=self.coord_frame,
            nside=self.nside,
        )
        map_data = HpMask(self.map_path, coordframe=self.coord_frame)

        demeaned_values = values - map_data.mask[cat_pix]
        catalog.add_column(demeaned_values, name="demeaned_k")

        return catalog

    def run(self):
        catalog = self.demean_catalog()
        catalog.write(self.output_path, format="fits", overwrite=True)
        return self.output_path
