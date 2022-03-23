"""
Another COS Dark Correction (ACDC)

1. Superdarks must be created. (Let's assume this is already done)
2. Bin superdarks and predict the dark level for each input corrtag.
3. Use the predicted dark level to apply a custom dark correction.
4. Calibrate the custom-corrected corrtags.
5. Coadd the calibrated x1ds to increase SNR.
"""

from collections import defaultdict
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
    def __init__(self, indir, outdir, binned=False):
        self.indir = indir
        self.outdir = outdir
        self.binned = binned 
        now = datetime.datetime.now()
        self.cal_outdir = os.path.join(outdir, f"cal_{now.strftime('%d%b%Y')}")
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        corrtags = glob.glob(os.path.join(indir, "*corrtag*fits"))
        self.corr_dict = self.sort_corrtags(corrtags)
        self.dark_dict = self.sort_superdarks()

    
    def sort_corrtags(self, corrtags):
        corr_dict = defaultdict(list)
        for item in corrtags:
            file_segment = fits.getval(item, "segment")
            file_hv = fits.getval(item, f"HVLEVEL{file_segment[-1]}", 1)
            corr_dict[f"{file_segment}_{file_hv}"].append(item)
        return corr_dict

    
    def sort_superdarks(self):
        darks0 = glob.glob(os.path.join(LOCAL_DARKDIR, "superdark*.asdf"))
        if self.binned is False:
            darks = [x for x in darks0 if "phabinned" not in x]
        else:
            darks = [x for x in darks0 if "phabinned" in x]
        dark_dict = defaultdict(dict)
        for dark in darks:
            darkfile = os.path.basename(dark)
            sp = darkfile.split("_")
            segment = sp[1]
            hv = sp[2]
            period = sp[-1].split(".")[0]
            dark_dict[f"{segment}_{hv}"][period] = dark
        return dark_dict


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
        for seg_hv in self.corr_dict:
            corrtags = self.corr_dict[seg_hv]
            lo_darkname = self.dark_dict[seg_hv]["quiescent"]
            hi_darkname = self.dark_dict[seg_hv]["active"]
            predict_dark(corrtags, lo_darkname, hi_darkname, 
                         outdir=self.outdir, binned=self.binned)
            subtract_dark(corrtags, self.outdir, outdir=self.outdir)
        self.custom_corrtags = glob.glob(os.path.join(self.outdir, "corrected*corrtag*fits"))

    
    def calibrate_corrtags(self):
        calibrated = []
        for item in self.custom_corrtags:
            if "_corrtag_a" in item:
                other = item.replace("_corrtag_a", "_corrtag_b")
            else:
                other = item.replace("_corrtag_b", "_corrtag_a")
            if other in calibrated:
                continue
            with fits.open(item, mode="update") as hdulist:
                hdr0 = hdulist[0].header
                hdr0.set("xtrctalg", "BOXCAR")
                hdr0.set("backcorr", "omit")
                hdr0.set("trcecorr", "omit")
                hdr0.set("algncorr", "omit")
            calcos.calcos(item, outdir=self.cal_outdir)
            calibrated.append(item)

    
#-----------------------------------------------------------------------------#

def run_acdc(indir, outdir, lo_darkname=None, hi_darkname=None, binned=False, 
             hv=None, segment=None):
    A = Acdc(indir, outdir, binned)
    #A = Acdc(indir, outdir, lo_darkname, hi_darkname, binned, hv, segment)
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
