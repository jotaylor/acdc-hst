import argparse
from acdc.dark_correction import Acdc


def run_acdc(indir, darkcorr_outdir, lo_darkname=None, hi_darkname=None, binned=False,
             hv=None, segment=None):
#TODO- allow for specific superdarks, hv, and segment
#TODO allow for specification of x1d_outdir
    """Wrapper script to perform custom dakr correction on an input directory.
    
    Args:
        indir (str): Input directory that houses corrtags to correct.
        darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
        lo_darkname (str): NOT YET IMPLEMENTED (Optional) Specific quiescent superdark to use. 
        hi_darkname (str): NOT YET IMPLEMENTED (Optional) Specific active superdark to use.
        binned (Bool): (Optional) If True, pre-binned superdarks will be used.
        hv (str): NOT YET IMPLEMENTED (Optional) Process corrtags of this HV only.
        segment (str): NOT YET IMPLEMENTED (Optional) Process corrtags of this segment only.
    """
    A = Acdc(indir=indir, darkcorr_outdir=darkcorr_outdir, binned=binned)
    #A = Acdc(indir, darkcorr_outdir, lo_darkname, hi_darkname, binned, hv, segment)
    A.custom_dark_correction()
    A.calibrate_corrtags()


def acdc_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",
                        help="Path to science corrtags")
    parser.add_argument("-o", "--darkcorr_outdir",
                        help="Name of output directory")
#    parser.add_argument("--lo", dest="lo_darkname", default=None,
#                        help="Name of low activity superdark")
#    parser.add_argument("--hi", dest="hi_darkname", default=None,
#                        help="Name of high activity superdark")
    parser.add_argument("--binned", default=False,
                        action="store_true",
                        help="Toggle to indicate that supplied superdarks are binned")
#    parser.add_argument("--hv", default=None,
#                        help="HV to filter corrtags by")
#    parser.add_argument("--segment", default=None,
#                        help="Segment to filter corrtags by")
    args = parser.parse_args()
#    run_acdc(args.indir, args.darkcorr_outdir, args.lo_darkname, args.hi_darkname, 
#             args.binned, args.hv, args.segment)
    run_acdc(indir=args.indir, darkcorr_outdir=args.darkcorr_outdir, binned=args.binned)


if __name__ == "__main__":
	acdc_parser()
