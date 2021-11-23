"""
Another COS Dark Correction (ACDC)

1. Superdarks must be created. (Let's assume this is already done)
2. Bin superdarks and predict the dark level for each input corrtag.
3. Use the predicted dark level to apply a custom dark correction.
4. Calibrate the custom-corrected corrtags.
5. Coadd the calibrated x1ds to increase SNR.
"""

import argparse
import os
import glob
import datetime
from astropy.io import fits
import calcos

from predict_dark_level1 import predict_dark
from subtr_standard import subtract_dark

# This only works for STScI internal folks
LOCAL_DARKDIR = "/astro/sveash/cos_dark/final_superdarks"

class Acdc():
    def __init__(self, indir, outdir, lo_darkname=None, hi_darkname=None, 
                 binned=False, hv=None, segment=None):
        self.indir = indir
        self.outdir = outdir
        now = datetime.datetime.now()
        self.cal_outdir = os.path.join(outdir, f"cal_{now.strftime('%d%b%Y')}")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        # To do handle determining darks programmatically    
        self.lo_darkname = lo_darkname
        self.hi_darkname = hi_darkname
        self.binned = binned 
        self.hv = hv
        self.segment = segment
        corrtags = glob.glob(os.path.join(indir, "*corrtag*fits"))
        if segment is not None or hv is not None:
            self.orig_corrtags = self.filter_corrtags(corrtags)
        else:
            self.orig_corrtags = corrtags


    def filter_corrtags(self, corrtags):
        good_corrtags = []
        for item in corrtags:
            file_segment = fits.getval(item, "segment")
            file_hv = fits.getval(item, f"HVLEVEL{file_segment[-1]}", 1)
            if self.segment is not None and file_segment != self.segment.upper():
                print(f"File does not match required settings: {item}")
                continue
            if self.hv is not None and int(file_hv) != int(self.hv):
                print(f"File does not match required settings: {item}")
                continue
            good_corrtags.append(item)
        assert len(good_corrtags) != 0, \
            "No supplied corrtags matched the settings HV={self.hv}, segment={self.segment}"
        return good_corrtags

# TODO- move binning out of custom correction
    def custom_dark_correction(self):
        predict_dark(self.orig_corrtags, self.lo_darkname, self.hi_darkname, 
                     outdir=self.outdir, binned=self.binned)
        subtract_dark(self.orig_corrtags, self.outdir, outdir=self.outdir)
        self.custom_corrtags = glob.glob(os.path.join(self.outdir, "corrected*corrtag*fits"))

    
    def calibrate_corrtags(self):
        for item in self.custom_corrtags:
            calcos.calcos(item, outdir=self.cal_outdir)


    

#-----------------------------------------------------------------------------#

def run_acdc(indir, outdir, lo_darkname=None, hi_darkname=None, binned=False, 
             hv=None, segment=None):
    A = Acdc(indir, outdir, lo_darkname, hi_darkname, binned, hv, segment)
    A.custom_dark_correction()
    A.calibrate_corrtags()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--indir",
                        help="Path to science corrtags")
    parser.add_argument("-o", "--outdir",
                        help="Name of output directory")
    parser.add_argument("--lo", dest="lo_darkname", default=None,
                        help="Name of low activity superdark")
    parser.add_argument("--hi", dest="hi_darkname", default=None,
                        help="Name of high activity superdark")
    parser.add_argument("--binned", default=False,
                        action="store_true",
                        help="Toggle to indicate that supplied superdarks are binned")
    parser.add_argument("--hv", default=None,
                        help="HV to filter corrtags by")
    parser.add_argument("--segment", default=None,
                        help="Segment to filter corrtags by")
    args = parser.parse_args()
    run_acdc(args.indir, args.outdir, args.lo_darkname, args.hi_darkname, 
             args.binned, args.hv, args.segment)
