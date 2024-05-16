import argparse
from acdc.dark_correction import Acdc


def run_acdc(indir, darkcorr_outdir, x1d_outdir=None,
             superdark_dir=None, binned=True, hv=None, 
             segment=None, overwrite=False, exclude_airglow=False,
             calibrate=True):
    """Wrapper script to perform custom dakr correction on an input directory.
    
    Args:
        indir (str): Input directory that houses corrtags to correct.
        darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
        x1d_outdir (str): Custom dark-corrected x1ds will be written here.
        binned (Bool): (Optional) If True, pre-binned superdarks will be used.
        hv (str): (Optional) Process corrtags of this HV only.
        segment (str): (Optional) Process corrtags of this segment only.
    """
    A = Acdc(indir=indir, darkcorr_outdir=darkcorr_outdir, x1d_outdir=x1d_outdir,
             superdark_dir=superdark_dir,
             binned=binned, segment=segment, hv=hv, overwrite=overwrite,
             exclude_airglow=exclude_airglow, calibrate=calibrate)
    A.custom_dark_correction()
    A.calibrate_corrtags()
    print(f"\nFINISHED: Products in {A.x1d_outdir}")

def acdc_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",
                        help="Path to science corrtags")
    parser.add_argument("-o", "--darkcorr_outdir",
                        help="Name of directory to write model superdarks and custom corrtags to")
    parser.add_argument("--x1d_outdir", default=None,
                        help="Name of directory to write 1D spectra to")
    parser.add_argument("-d", "--darkdir", default=None,
                        help="Directory that houses superdarks to use")
    parser.add_argument("--unbinned", dest="binned", default=True,
                        action="store_false",
                        help="Toggle to indicate that supplied superdarks are binned")
    parser.add_argument("-c", "--clobber", default=False,
                        action="store_true",
                        help="Toggle to overwrite any existing products")
    parser.add_argument("--hv", default=None,
                        help="HV to filter corrtags by")
    parser.add_argument("--segment", default=None,
                        help="Segment to filter corrtags by")
    parser.add_argument("--airglow", default=False, dest="exclude_airglow",
                        action="store_true",
                        help="Toggle to exclude LyA and OI airglow from science data scaling") 
    parser.add_argument("--nocal", default=True, dest="calibrate",
                        action="store_false",
                        help="Toggle to turn off CalCOS calibration") 
    args = parser.parse_args()
    run_acdc(args.indir, args.darkcorr_outdir, args.x1d_outdir, 
             args.darkdir, args.binned, args.hv, 
             args.segment, args.clobber, args.exclude_airglow, args.calibrate)


if __name__ == "__main__":
	acdc_parser()
