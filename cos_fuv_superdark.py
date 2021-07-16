import datetime
import argparse
import os
import glob

import asdf
from astropy.io import fits
from astropy.table import Table
from calcos import ccos

import numpy as np
import pandas as pd

from query_cos_dark import files_by_mjd

class Superdark():
    def __init__(self, hv, segment, mjdstarts, mjdends, dayint=100, bin_x=1,
                 bin_y=1, bin_pha=1, phastart=1, phaend=31, 
                 gsagtab="41g2040ol_gsag.fits", region="inner", outfile=None):

        if segment == "FUVA":
            if region == "inner":
                xstart, xend, ystart, yend = 1260, 15119, 375, 660
        elif segment == "FUVB":
            if region == "inner":
                xstart, xend, ystart, yend = 1000, 14990, 405, 740
        else:
            raise Exception(f"Invalid segment specified: {segment}")
        self.xstart, self.ystart = xstart, ystart
        self.xend = (xend // bin_x) * bin_x
        self.yend = (yend // bin_y) * bin_y
        self.bin_xstart = self.xstart // bin_x
        self.bin_ystart = self.ystart // bin_y
        self.bin_xend = self.xend // bin_x
        self.bin_yend = self.yend // bin_y

        self.segment = segment
        self.hv = hv
        self.region = region
        self.mjdstarts = mjdstarts
        self.mjdends = mjdends
        self.ndays = [mjdends[i] - mjdstarts[i] for i in range(len(mjdstarts))]
        self.dayint = dayint
        self.bin_x = bin_x
        self.bin_y = bin_y
        if bin_pha == "all":
            bin_pha = phaend - phastart - 1
        self.bin_pha = bin_pha
        self.gsagtab = gsagtab
        self.superdark_unbinned = np.zeros(1024*16384).reshape((1024, 16384))
        self.outfile = outfile

        self.get_gsag_holes()


    def create_superdark(self):
        notfilled = True
        pha_images = {}
        for k in range(len(self.mjdstarts)):
            start = self.mjdstarts[k]
            total_days = 0
            total_exptime = 0
            total_files = 0
            pha_range = np.arange(1, 31, self.bin_pha)
            if pha_range[-1] != 30:
                pha_range = np.concatenate( (pha_range, np.array([30])) )

            runstart = datetime.datetime.now()
            print("Start time: {}".format(runstart))
            while notfilled is True:
                total_days += self.dayint
                if total_days > self.ndays[k]:
                    total_days = self.ndays[k]
                end = start + self.dayint
                if end > self.mjdends[k]:
                    end = self.mjdends[k]
                notfilled = False
                print(f"Using darks from MJD {start:,}-{end:,}; {total_days}/{self.ndays[k]} days")
                gsag_df = self.gsag_holes.loc[(self.gsag_holes["DATE"] > start) & (self.gsag_holes["DATE"] < end)]
                print("   Querying...")
                dark_df = files_by_mjd(start, end, segment=self.segment, hv=self.hv)
                print("   Query done")
                darks = dark_df["fileloc"].values
                if len(darks) == 0:
                    notfilled = True
                    continue
                total_exptime += sum([fits.getval(x, "exptime", 1) for x in darks])
                total_files += len(darks)
                for i in range(len(pha_range)-1):
                    print("   Binning corrtags...")
                    sum_image = self.bin_corrtags(darks, phastart=pha_range[i], phaend=pha_range[i+1])
                    print("   Binning done")
                    for j in range(len(gsag_df)):
                        sum_image[gsag_df.iloc[j]["Y0"]:gsag_df.iloc[j]["Y1"], gsag_df.iloc[j]["X0"]:gsag_df.iloc[j]["X1"]] = 99999
                    tmp = sum_image.reshape(1024 // self.bin_y, self.bin_y,
                                            16384 // self.bin_x, self.bin_x)
                    binned = tmp.sum(axis=3).sum(axis=1)
                    binned_inner = binned[self.bin_ystart:self.bin_yend, self.bin_xstart:self.bin_xend]
                    zeros = np.where(binned_inner == 0)
                    if len(zeros[0]) != 0:
                        notfilled = True
                    key = f"pha{pha_range[i]}-{pha_range[i+1]}"
                    if key in pha_images:
                        pha_images[key] += binned_inner
                    else:
                        pha_images[key] = binned_inner
                start = end
                if total_days >= self.ndays[k]:
                    print("\nWarning: Not every pixel had events at every PHA")
                    notfilled = False
        # Undo 99999 put in for gainsag holes
        for i in range(len(pha_range)-1):
            key = f"pha{pha_range[i]}-{pha_range[i+1]}"
            inds = np.where(pha_images[key] == 99999)
            pha_images[key][inds] = 0
        pha_images["xstart"] = self.xstart
        pha_images["ystart"] = self.ystart
        pha_images["phastart"] = pha_range[0]
        pha_images["xend"] = self.xend
        pha_images["yend"] = self.yend
        pha_images["phaend"] = pha_range[-1]
        pha_images["bin_x"] = self.bin_x
        pha_images["bin_y"] = self.bin_y
        pha_images["bin_pha"] = self.bin_pha
        pha_images["mjdstarts"] = self.mjdstarts
        pha_images["mjdends"] = self.mjdends
        pha_images["segment"] = self.segment
        pha_images["hv"] = self.hv
        pha_images["total_exptime"] = total_exptime
        pha_images["total_files"] = total_files
        af = asdf.AsdfFile(pha_images)
        today = runstart.strftime("%d%b%y-%H%M")
        if self.outfile is None:
            self.outfile = f"superdark_{self.segment}_{self.hv}_{today}.asdf"
        af.write_to(self.outfile)
        print(f"Wrote {self.outfile}")
        runend = datetime.datetime.now()
        print("End time: {}".format(runend))
        print("Total time: {}".format(runend-runstart))


    def get_gsag_holes(self):
        extfound = False
        extnum = 1
        hvkey = f"HVLEVEL{self.segment[-1]}"
        while extfound is False:
            try:
                gsag_hv = fits.getval(self.gsagtab, hvkey, extnum)
                gsag_seg = fits.getval(self.gsagtab, "segment", extnum)
            except KeyError:
                extnum += 1
                continue
            if gsag_hv == self.hv and gsag_seg == self.segment:
                extfound = True
            else:
                extnum += 1
        gsag = Table.read(self.gsagtab, format="fits", hdu=extnum)
        gsag_df = gsag.to_pandas()
        gsag_df = gsag_df.rename(columns={"LX": "X0", "LY": "Y0"})
        gsag_df["X1"] = gsag_df["X0"] + gsag_df["DX"]
        gsag_df["Y1"] = gsag_df["Y0"] + gsag_df["DY"]
        gsag_df["BIN_X0"] = gsag_df["X0"]//self.bin_x
        gsag_df["BIN_X1"] = gsag_df["X1"]//self.bin_x
        gsag_df["BIN_Y0"] = gsag_df["Y0"]//self.bin_y
        gsag_df["BIN_Y1"] = gsag_df["Y1"]//self.bin_y
        gsag_df = gsag_df.astype("int32")

        self.gsag_holes = gsag_df

    def bin_corrtags(self, corrtag_list, phastart, phaend, xtype='XCORR', ytype='YCORR', sdqflags=0):
        """Bin one or more corrtag event lists into a 2D image.

        Modifed from /user/jotaylor/git/jotools/utils.py

        Parameters
        ----------
        corrtag_list : str
            list of datasets to combine into an image
        xtype : str, optional
            X-coordinate type
        ytype : str, optional
            Y-coordinate type

        Returns
        -------
        image : np.ndarray
            2D image of the corrtag
        """

        final_image = np.zeros((1024, 16384)).astype(np.float32)

        for filename in corrtag_list:
            image = np.zeros((1024, 16384)).astype(np.float32)
            events = fits.getdata(filename, 1)
            inds = np.where((events["pha"] >= phastart) & (events["pha"] < phaend))
            events = events[inds]
            if not len(inds[0]):
                print(f"No counts for PHA={pha}")
                return image

            # Call for this is x_values, y_values, image to bin to, offset in x
            # ccos.binevents(x, y, array, x_offset, dq, sdqflags, epsilon)
            ccos.binevents(events[xtype].astype(np.float32),
                           events[ytype].astype(np.float32),
                           image,
                           0,
                           events['dq'],
                           sdqflags)

            final_image += image

        return final_image


