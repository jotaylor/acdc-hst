"""
Perform a custom dark correction on COS/FUV data.
For each input science exposure, a quiescent and an active superdark 
of the appropriate segment+HV combination are used to determine the dark rate 
of the science exposure. A custom dark correction is then applied using
this model dark rate, creating custom corrtag products. These corrtags are
calibrated using CalCOS to create custom x1d files. x1d files should be coadded
offline in order to increase SNR.
"""

import sys
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
                 superdark_dir=None, segment=None, hv=None, overwrite=False,
                 exclude_lya=False, lo_darkname=None, hi_darkname=None):
        """
        Args:
            indir (str): Input directory that houses corrtags to correct.
            darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
            x1d_outdir (str): Custom dark-corrected x1ds will be written here.
            binned (Bool): If True, pre-binned superdarks will be used.
            superdark_dir (str): Location of superdarks. 
        """

        self.overwrite = overwrite
        self.indir = indir
        self.exclude_lya = exclude_lya
        if superdark_dir is None:
            try:
                superdark_dir = os.environ["ACDC_SUPERDARKS"]
            except KeyError as e:
                print("ERROR: You must define the $ACDC_SUPERDARKS environment variable- this is where all superdarks are located")
                print("Exiting")
                sys.exit()
        self.superdark_dir = superdark_dir

        self.darkcorr_outdir = darkcorr_outdir
        self.binned = binned
        if segment is not None:
            self.segment = segment.upper()
        else:
            self.segment = segment
        self.hv = hv 
        now = datetime.datetime.now()
        self.x1d_outdir = os.path.join(darkcorr_outdir, f"cal_{now.strftime('%d%b%Y')}")
        if not os.path.exists(darkcorr_outdir):
            os.makedirs(darkcorr_outdir)
        corrtags = glob.glob(os.path.join(indir, "*corrtag*fits"))
        self.corr_dict = self.sort_corrtags(corrtags)

        if set([lo_darkname, hi_darkname]) == {None}:
            supplied_darks = None
        if None in [lo_darkname, hi_darkname] and len(set([lo_darkname, hi_darkname])) != 1:
            print(f"WARNING: Both an active and quiescent must be provided- using default superdark library")
            supplied_darks = None
        elif len(set([lo_darkname, hi_darkname])) == 2:
            supplied_darks = [lo_darkname, hi_darkname]
        self.dark_dict = self.sort_superdarks(supplied_darks)

    
    def sort_corrtags(self, corrtags):
        """
        Sort corrtags into a dictionary. If specified, filter corrtags based 
        on segment and HV setting.
        
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
            if self.segment is not None:
                if file_segment != self.segment:
                    continue
            if self.hv is not None:
                if file_hv != self.hv:
                    continue
            corr_dict[f"{file_segment}_{file_hv}"].append(item)
        return corr_dict

    
    def sort_superdarks(self, supplied_darks=None):
        """Sort superdarks into a dictionary based on segment and HV setting.
        
        Returns:
            dark_dict (dict): Nested dictionary where each key is the segment+HV setting
                and each value is a dictionary where each key is the type of superdark 
                (either 'active' or 'quiescent') and each value is the superdark name.
        """

#TODO - handle multiple superdarks per activity period
        if supplied_darks is None:
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
                         outdir=self.darkcorr_outdir, binned=self.binned,
                         overwrite=self.overwrite, segment=self.segment,
                         hv=self.hv, exclude_lya=self.exclude_lya)
            subtract_dark(corrtags, self.darkcorr_outdir, outdir=self.darkcorr_outdir, 
                          overwrite=self.overwrite)
        self.custom_corrtags = glob.glob(os.path.join(self.darkcorr_outdir, "corrected*corrtag*fits"))

    
    def calibrate_corrtags(self):
        """Calibrate custom dark-corrected corrtags to product x1ds.

        Calibrate each custom dark-corrected corrtag using the BOXCAR extraction
        method in CalCOS. The TWOZONE extraction method does not work with
        BACKCORR=OMIT.
        """

        for item in self.custom_corrtags:
            with fits.open(item, mode="update") as hdulist:
                hdr0 = hdulist[0].header
                hdr0.set("xtrctalg", "BOXCAR")
                hdr0.set("backcorr", "OMIT")
                hdr0.set("trcecorr", "OMIT")
                hdr0.set("algncorr", "OMIT")

        calibrated = []
        for item in self.custom_corrtags:
            if "_corrtag_a" in item:
                other = item.replace("_corrtag_a", "_corrtag_b")
            else:
                other = item.replace("_corrtag_b", "_corrtag_a")
            if other in calibrated:
                continue
            spl = os.path.basename(item).split("_")
            wildcard = spl[0] + "_" + spl[1] + "*"
            wildcard_products = glob.glob(os.path.join(self.x1d_outdir, wildcard))
            if len(wildcard_products) >= 1 and self.overwrite is True:
                for item in wildcard_products:
                    os.remove(item)
            elif len(wildcard_products) >= 1 and self.overwrite is False:
                print(f"WARNING: Products already exist and overwrite is False, skipping...")
                continue

            calcos.calcos(item, outdir=self.x1d_outdir, verbosity=0)
            calibrated.append(item)
            calibrated.append(other)

