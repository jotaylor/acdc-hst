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

from .predict_dark_level1 import predict_dark
from .subtr_standard import subtract_dark

# This only works for STScI internal folks
LOCAL_DARKDIR = "/astro/sveash/cos_dark/final_superdarks"

class Acdc():
    """Perform a custom dark correction on COS/FUV data.

    Attributes:
        indir (str): Input directory that houses corrtags to correct.
        darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
        x1d_outdir (str): Custom dark-corrected x1ds will be written here.
        binned (Bool): If True, pre-binned superdarks will be used.
        corr_dict (dict): Dictionary where each key is the segment+HV setting 
            and each value is a list of all corrtags with that setting.
        custom_corrtags (list): List of corrtags which have had the custom
            dark correction applied.
        dark_dict (dict): Nested dictionary where each key is the segment+HV setting
            and each value is a dictionary where each key is the type of superdark 
            (either 'active' or 'quiescent') and each value is the superdark name.
    """
    
    def __init__(self, indir, darkcorr_outdir, x1d_outdir=None, binned=False, 
                 superdark_dir=LOCAL_DARKDIR):
        """
        Args:
            indir (str): Input directory that houses corrtags to correct.
            darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
            x1d_outdir (str): Custom dark-corrected x1ds will be written here.
            binned (Bool): If True, pre-binned superdarks will be used.
            superdark_dir (str): Location of superdarks. 
        """

        self.indir = indir
        self.superdark_dir = superdark_dir
        self.darkcorr_outdir = darkcorr_outdir
        self.binned = binned 
        now = datetime.datetime.now()
        self.x1d_outdir = os.path.join(darkcorr_outdir, f"cal_{now.strftime('%d%b%Y')}")
        if not os.path.exists(darkcorr_outdir):
            os.makedirs(darkcorr_outdir)
        corrtags = glob.glob(os.path.join(indir, "*corrtag*fits"))
        self.corr_dict = self.sort_corrtags(corrtags)
        self.dark_dict = self.sort_superdarks()

    
    def sort_corrtags(self, corrtags):
        """Sort corrtags into a dictionary based on segment and HV setting.
        
        Args:
            corrtags (array-like): All corrtags to be corrected.
        
        Returns:
        corr_dict (dict): Dictionary where each key is the segment+HV setting 
            and each value is a list of all corrtags with that setting.
        """
        
        corr_dict = defaultdict(list)
        for item in corrtags:
            file_segment = fits.getval(item, "segment")
            file_hv = fits.getval(item, f"HVLEVEL{file_segment[-1]}", 1)
            corr_dict[f"{file_segment}_{file_hv}"].append(item)
        return corr_dict

    
    def sort_superdarks(self):
        """Sort superdarks into a dictionary based on segment and HV setting.
        
        Returns:
        dark_dict (dict): Nested dictionary where each key is the segment+HV setting
            and each value is a dictionary where each key is the type of superdark 
            (either 'active' or 'quiescent') and each value is the superdark name.
        """

#TODO - handle multiple superdarks per activity period
        darks0 = glob.glob(os.path.join(self.superdark_dir, "superdark*.asdf"))
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
        """Filter corrtags based on HV and segment requirements.

        Args:
            corrtags (array-like): All corrtags to be corrected.

        Returns:
            good_corrtags (array-like): All corrtags to be corrected,
                filtered by the specified HV and segment requirements.
        """
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
        """Perform the custom dark correction.

        For each corrtag, use the superdarks to predict the dark rate across 
        the entire detector. Then use this predicted dark rate to perform a 
        custom dark correction to each corrtag.
        """
        
        for seg_hv in self.corr_dict:
            corrtags = self.corr_dict[seg_hv]
            lo_darkname = self.dark_dict[seg_hv]["quiescent"]
            hi_darkname = self.dark_dict[seg_hv]["active"]
            predict_dark(corrtags, lo_darkname, hi_darkname, 
                         outdir=self.darkcorr_outdir, binned=self.binned)
            subtract_dark(corrtags, self.darkcorr_outdir, outdir=self.darkcorr_outdir)
        self.custom_corrtags = glob.glob(os.path.join(self.darkcorr_outdir, "corrected*corrtag*fits"))

    
    def calibrate_corrtags(self):
        """Calibrate custom dark-corrected corrtags to product x1ds.

        Calibrate each custom dark-corrected corrtag using the BOXCAR extraction
        method in CalCOS. The TWOZONE extraction method does not work with
        BACKCORR=OMIT.
        """

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
            calcos.calcos(item, outdir=self.x1d_outdir)
            calibrated.append(item)

    
#-----------------------------------------------------------------------------#

def run_acdc(indir, darkcorr_outdir, lo_darkname=None, hi_darkname=None, binned=False, 
             hv=None, segment=None):
#TODO- allow for specific superdarks, hv, and segment
#TODO allow for specification of x1d_outdir
    """Wrapper script to perform custom dakr correction on an input directory.
    
    Args:
        indir (str): Input directory that houses corrtags to correct.
        darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
        binned (Bool): If True, pre-binned superdarks will be used.
    """
    A = Acdc(indir, outdir, binned)
    #A = Acdc(indir, outdir, lo_darkname, hi_darkname, binned, hv, segment)
    A.custom_dark_correction()
    A.calibrate_corrtags()

if __name__ == "__main__":
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
    run_acdc(args.indir, args.darkcorr_outdir, args.binned)
