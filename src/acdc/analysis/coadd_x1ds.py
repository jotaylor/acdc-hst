import os
import glob
import argparse
import datetime

from ullyses.generic_coadd_wrapper import coadd_files

def coadd_cos_spectra(infiles, outdir=None, outfile=None, clobber=False):
    """Coadd multiple 1-D COS spectra of the same grating.
    Args:
        grating (str): COS grating name.
        indir (str): Path that holds all data to be coadded. All data in this
            directory that match grating will be coadded.
    Returns:
        None
    """

    if outdir is None:
        now = datetime.datetime.now()
        outdir = f"{now.strftime('%d%b%Y_%H%M%S')}_coadd"  
    coadd_files(infiles, outdir, outfile, clobber)


def find_files(indir):
    allfiles = glob.glob(os.path.join(indir, "*x1d.fits"))
    return allfiles

def coadd_parser():
    """
    Copied from ullyses.ullyses_coadd_abut_wrapper
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",
                        default="./",
                        help="Directory(ies) with data to combine")
    parser.add_argument("-o", "--outdir", default=None,
                        help="Directory for output HLSPs")
    parser.add_argument("--outfile", default=None,
                        help="Name of output coadded file")
    parser.add_argument("-c", "--clobber", default=False,
                        action="store_true",
                        help="If True, overwrite existing products")
    args = parser.parse_args()

    infiles = find_files(args.indir)
    coadd_cos_spectra(infiles, args.outdir, args.outfile, args.clobber)


if __name__ == "__main__":
    coadd_parser()
