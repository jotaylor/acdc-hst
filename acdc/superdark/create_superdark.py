import numpy as np
import argparse

from acdc.superdark.cos_fuv_superdark import Superdark

def all_hvs():
    S = Superdark(hv=163, segment= "FUVA", mjdstarts=[58000], mjdends=[58220], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_163_58000_58220_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=163, segment= "FUVA", mjdstarts=[58300], mjdends=[58600], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_163_58300_58600_active.asdf")
    S.create_superdark()
    S = Superdark(hv=167, segment= "FUVA", mjdstarts=[57500], mjdends=[57800], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_167_57500_57800_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=167, segment= "FUVA", mjdstarts=[56520], mjdends=[56820], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_167_56520_56820_active.asdf")
    S.create_superdark()
    S = Superdark(hv=169, segment= "FUVA", mjdstarts=[55100], mjdends=[55400], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_169_55100_55400_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=169, segment= "FUVA", mjdstarts=[55700], mjdends=[56010], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_169_55700_56010_active.asdf")
    S.create_superdark()
    S = Superdark(hv=171, segment= "FUVA", mjdstarts=[57050], mjdends=[57300], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_171_57050_57300_active.asdf")
    S.create_superdark()
    S = Superdark(hv=173, segment= "FUVA", mjdstarts=[56900], mjdends=[57100], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVA_173_56900_57100_active.asdf")
    S.create_superdark()
    S = Superdark(hv=163, segment= "FUVB", mjdstarts=[56250], mjdends=[56420], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_163_56250_56420_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=163, segment= "FUVB", mjdstarts=[58000], mjdends=[58300], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_163_58000_58300_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=163, segment= "FUVB", mjdstarts=[57050], mjdends=[57250], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_163_57050_57250_active.asdf")
    S.create_superdark()
    S = Superdark(hv=167, segment= "FUVB", mjdstarts=[55100], mjdends=[55400], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_167_55100_55400_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=167, segment= "FUVB", mjdstarts=[57000], mjdends=[57300], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_167_57000_57300_active.asdf")
    S.create_superdark()
    S = Superdark(hv=169, segment= "FUVB", mjdstarts=[57400], mjdends=[57700], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_169_57400_57700_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=169, segment= "FUVB", mjdstarts=[56550], mjdends=[56800], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_169_56550_56800_active.asdf")
    S.create_superdark()
    S = Superdark(hv=175, segment= "FUVB", mjdstarts=[57650], mjdends=[58050], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_175_57650_58050_quiescent.asdf")
    S.create_superdark()
    S = Superdark(hv=175, segment= "FUVB", mjdstarts=[56850], mjdends=[57100], bin_pha=1, phastart=1, phaend=31, 
                  outfile="superdark_FUVB_175_56850_57100_active.asdf")
    S.create_superdark()

if __name__ == "__main__":
#    parser = argparse.ArgumentParser()
#    parser.add_argument("--hv", type=int, help="HV setting of interest")
#    parser.add_argument("--segment", help="Segment of interest")
#    parser.add_argument("--mjdstarts", nargs="+", help="MJD start dates")
#    parser.add_argument("--mjdends", nargs="+", help="MJD end dates") 
#    parser.add_argument("--phabin", default=1,
#                        help="Size of PHA binning")
#    parser.add_argument("--phastart", default=1, type=int,
#                        help="Minimum PHA to consider")
#    parser.add_argument("--phaend", default=31, type=int,
#                        help="Maximum PHA to consider")
#    parser.add_argument("--pha_bins", nargs="*", default=None,
#                        help="Explicit PHA bins") 
#    parser.add_argument("-o", "--outfile", default=None,
#                        help="Name of output superdark")
#    args = parser.parse_args()
#
#    mjdstarts = [int(x) for x in args.mjdstarts]
#    mjdends = [int(x) for x in args.mjdends]
#    if args.pha_bins is not None:
#        pha_bins = np.array([int(x) for x in args.pha_bins])
#    S = Superdark(hv=args.hv, segment=args.segment, mjdstarts=mjdstarts, 
#                  mjdends=mjdends, bin_pha=args.phabin, phastart=args.phastart, 
#                  phaend=args.phaend, pha_bins=pha_bins, outfile=args.outfile)
#    S.create_superdark()
    all_hvs()
