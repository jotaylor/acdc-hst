import numpy as np
import argparse

from cos_fuv_superdark import Superdark

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--hv", type=int, help="HV setting of interest")
    parser.add_argument("--segment", help="Segment of interest")
    parser.add_argument("--mjdstarts", nargs="+", help="MJD start dates")
    parser.add_argument("--mjdends", nargs="+", help="MJD end dates") 
    parser.add_argument("--phabin", default=1,
                        help="Size of PHA binning")
    parser.add_argument("--phastart", default=1, type=int,
                        help="Minimum PHA to consider")
    parser.add_argument("--phaend", default=31, type=int,
                        help="Maximum PHA to consider")
    parser.add_argument("--pha_bins", nargs="*", default=None,
                        help="Explicit PHA bins") 
    parser.add_argument("-o", "--outfile", default=None,
                        help="Name of output superdark")
    args = parser.parse_args()

    mjdstarts = [int(x) for x in args.mjdstarts]
    mjdends = [int(x) for x in args.mjdends]
    if args.pha_bins is not None:
        pha_bins = np.array([int(x) for x in args.pha_bins])
    S = Superdark(hv=args.hv, segment=args.segment, mjdstarts=mjdstarts, 
                  mjdends=mjdends, bin_pha=args.phabin, phastart=args.phastart, 
                  phaend=args.phaend, pha_bins=pha_bins, outfile=args.outfile)
    S.create_superdark()

