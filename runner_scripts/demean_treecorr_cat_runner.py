import argparse

from src.demean_treecorr_cat import DemeanTreecorrCat


DEFAULT_CAT_DIR = "/n23data1/mccleary/dustyhalos/dusthalos/output/sdss_csfd"
DEFAULT_CAT_NAME = "dustcorrel_bgr_photoz2_treecorrcat.fits"
DEFAULT_MAP_FILE = "masked_sdss_csfd_av_map_rand.fits"
DEFAULT_OUTPUT_NAME = "demeaned_bgr_photoz2_treecorrcat.fits"
DEFAULT_NSIDE = 256


def parse_args():
    parser = argparse.ArgumentParser(
        description="De-mean a TreeCorr catalog using a HEALPix map."
    )
    parser.add_argument(
        "--catdir",
        type=str,
        default=DEFAULT_CAT_DIR,
        help="Directory containing the input catalog and HEALPix mask.",
    )
    parser.add_argument(
        "--catname",
        type=str,
        default=DEFAULT_CAT_NAME,
        help="Filename of the input TreeCorr catalog FITS file.",
    )
    parser.add_argument(
        "--mapfile",
        type=str,
        default=DEFAULT_MAP_FILE,
        help="Filename of the HEALPix mask FITS file in cat-dir.",
    )
    parser.add_argument(
        "--outname",
        type=str,
        default=DEFAULT_OUTPUT_NAME,
        help="Filename to save the demeaned catalog as.",
    )
    parser.add_argument(
        "--nside",
        type=int,
        default=DEFAULT_NSIDE,
        help="HEALPix NSIDE used to match catalog rows to map pixels.",
    )
    parser.add_argument(
        "--ra-col",
        type=str,
        default="ra",
        help="Catalog column containing right ascension.",
    )
    parser.add_argument(
        "--dec-col",
        type=str,
        default="dec",
        help="Catalog column containing declination.",
    )
    parser.add_argument(
        "--value-col",
        type=str,
        default="k",
        help="Catalog column containing the quantity to be demeaned.",
    )
    parser.add_argument(
        "--coordframe",
        type=str,
        default="icrs",
        help="Coordinate frame used for the HEALPix conversion, e.g. icrs or galactic.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    demeaner = DemeanTreecorrCat(
        cat_dir=args.cat_dir,
        cat_name=args.cat_name,
        map_file=args.map_file,
        output_name=args.output_name,
        nside=args.nside,
        ra_col=args.ra_col,
        dec_col=args.dec_col,
        value_col=args.value_col,
        coord_frame=args.coord_frame,
    )
    return demeaner.run()


if __name__ == "__main__":
    main()
