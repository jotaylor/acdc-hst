"""
Another COS Dark Correction (ACDC)
Perform a custom dark correction on COS/FUV data.
For each input science exposure, a quiescent and an active superdark 
of the appropriate segment+HV combination are used to determine the dark rate 
of the science exposure. A custom dark correction is then applied using
this model dark rate, creating custom corrtag products. These corrtags are
calibrated using CalCOS to create custom x1d files. x1d files should be coadded
offline in order to increase SNR.
"""

from collections import defaultdict
import argparse
import os
import glob
import datetime
from astropy.io import fits
import calcos

from acdc.predict_dark_level1 import predict_dark
from acdc.subtr_standard import subtract_dark


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
                 superdark_dir=None):
        """
        Args:
            indir (str): Input directory that houses corrtags to correct.
            darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
            x1d_outdir (str): Custom dark-corrected x1ds will be written here.
            binned (Bool): If True, pre-binned superdarks will be used.
            superdark_dir (str): Location of superdarks. 
        """

        self.indir = indir
        if superdark_dir is None:
            try:
                superdark_dir = os.environ["ACDC_SUPERDARKS"]
            except KeyError as e:
                print(e.message)
                print("You must define the $ACDC_SUPERDARKS environment variable- this is where all superdarks are located")
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
                hdr0.set("backcorr", "OMIT")
                hdr0.set("trcecorr", "OMIT")
                hdr0.set("algncorr", "OMIT")
            calcos.calcos(item, outdir=self.x1d_outdir)
            calibrated.append(item)

    

