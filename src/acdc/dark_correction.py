"""
Perform a custom dark correction on COS/FUV data.
For each input science exposure, a quiescent and an active superdark 
of the appropriate segment+HV combination are used to determine the dark rate 
of the science exposure. A custom dark correction is then applied using
this model dark rate, creating custom corrtag products. These corrtags are
calibrated using CalCOS to create custom x1d files. x1d files should be coadded
offline in order to increase SNR.
"""

import os
DIRNAME = os.path.dirname(__file__)
import sys
import gzip
import shutil
from collections import defaultdict
import argparse
import glob
import datetime
from astropy.io import fits
import asdf
import calcos

from acdc.predict_dark_level import predict_dark
from acdc.subtract_dark_level import subtract_dark


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
    
    def __init__(self, indir, darkcorr_outdir, x1d_outdir=None, binned=True, 
                 segment=None, hv=None, overwrite=False,
                 exclude_airglow=False, superdark_dir=None, 
                 superdarks=None, calibrate=True):
        """
        Args:
            indir (str): Input directory that houses corrtags to correct.
            darkcorr_outdir (str): Custom dark-corrected corrtags will be written here.
            x1d_outdir (str): Custom dark-corrected x1ds will be written here.
            binned (Bool): If True, pre-binned superdarks will be used.
            superdark_dir (str): Location of superdarks. 
        """

        self.calibrate = calibrate
        self.overwrite = overwrite
        self.indir = indir
        self.exclude_airglow = exclude_airglow
        if superdark_dir is None and superdarks is None:
            trydir = os.path.join(DIRNAME, "data/superdarks/epoch_superdarks")
            if len(glob.glob(os.path.join(trydir, "*superdark*asdf*"))) != 0:
                superdark_dir = trydir 
            else:
                try:
                    superdark_dir = os.environ["ACDC_SUPERDARKS"]
                except KeyError as e:
                    print("ERROR: You must supply the superdark directory or define the $ACDC_SUPERDARKS environment variable- this is where all superdarks are located")
                    print("Exiting")
                    sys.exit()
        self.superdark_dir = superdark_dir
        if superdarks is None:
            superdarks = glob.glob(os.path.join(self.superdark_dir, "*superdark*.asdf*"))
#        superdarks = self.unzip_superdarks(superdarks)
        
        self.darkcorr_outdir = darkcorr_outdir
        self.binned = binned
        if segment is not None:
            segment = segment.upper()
        self.segment = segment 
        if hv is not None:
            hv = int(hv)
        self.hv = hv 
        now = datetime.datetime.now()
        self.x1d_outdir = os.path.join(darkcorr_outdir, f"cal_{now.strftime('%d%b%Y')}")
        if not os.path.exists(darkcorr_outdir):
            os.makedirs(darkcorr_outdir)
        corrtags = glob.glob(os.path.join(indir, "*corrtag*fits"))
        self.corr_dict = self.sort_corrtags(corrtags)
        self.dark_dict = self.get_best_superdarks(superdarks)

    
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
                if int(file_hv) != self.hv:
                    continue
            corr_dict[f"{file_segment}_{file_hv}"].append(item)
        assert len(corr_dict) != 0, "No corrtags found that match segment and HV constraints"
        return corr_dict

    
    def get_best_superdarks(self, superdarks):
        """Sort superdarks into a dictionary based on segment and HV setting.
        
        Returns:
            dark_dict (dict): Dictionary where each key is the segment+HV setting
                and each value is a list of all applicable superdarks. 
        """
        if self.binned is False:
            darks = [x for x in superdarks if "binned" not in x]
        else:
            darks = [x for x in superdarks if "binned" in x]
        dark_dict = defaultdict(list)
        corr_segments = [x.split("_")[0] for x in self.corr_dict]
        corr_hvs = [x.split("_")[1] for x in self.corr_dict]
        for dark in darks:
            darkfile = os.path.basename(dark)
            sp = darkfile.split("_")
            segment = sp[1]
            hv = sp[2]
            for i in range(len(corr_hvs)):
                if corr_segments[i] == segment and corr_hvs[i] == hv:
                    dark_dict[f"{segment}_{hv}"].append(dark)
        
        assert len(dark_dict) != 0, "No matching superdarks found!!"
        if self.binned is False:
            print("Matching unbinned superdarks:")
        else:
            print("Matching binned superdarks:")
        for k,v in dark_dict.items():
            for sd in v:
                print(f"\t{sd}")

        return dark_dict

    
    def unzip_superdarks(self, superdarks):
        unzipped_superdarks = []
        for item in superdarks:
            if item.endswith(".gz"):
                unzipped = item.replace(".gz", "")
                with gzip.open(item, 'rb') as f_in:
                    with open(unzipped, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                print(f"Unzipped {item}")
                unzipped_superdarks.append(unzipped)
            else:
                unzipped_superdarks.append(item)
        return unzipped_superdarks
   

    def custom_dark_correction(self):
        """Perform the custom dark correction.

        For each corrtag, use the superdarks to predict the dark rate across 
        the entire detector. Then use this predicted dark rate to perform a 
        custom dark correction to each corrtag.
        """
        
        for seg_hv in self.corr_dict:
            corrtags = self.corr_dict[seg_hv]
            superdarks = self.dark_dict[seg_hv]
            predict_dark(corrtags, superdarks, 
                         outdir=self.darkcorr_outdir, binned=self.binned,
                         overwrite=self.overwrite, segment=self.segment,
                         hv=self.hv, exclude_airglow=self.exclude_airglow)
            subtract_dark(corrtags, self.darkcorr_outdir, outdir=self.darkcorr_outdir, 
                          overwrite=self.overwrite)
        self.custom_corrtags = glob.glob(os.path.join(self.darkcorr_outdir, "corrected*corrtag*fits"))

    
    def calibrate_corrtags(self):
        """Calibrate custom dark-corrected corrtags to product x1ds.

        Calibrate each custom dark-corrected corrtag using the BOXCAR extraction
        method in CalCOS. The TWOZONE extraction method does not work with
        BACKCORR=OMIT.
        """

        if self.calibrate is False:
            print("Calibration is set to skip, ending ACDC now")
            return

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
                for prod in wildcard_products:
                    os.remove(prod)
            elif len(wildcard_products) >= 1 and self.overwrite is False:
                print(f"WARNING: Products already exist and overwrite is False, skipping...")
                continue

            calcos.calcos(item, outdir=self.x1d_outdir, verbosity=0)
            calibrated.append(item)
            calibrated.append(other)

