import os
os.environ["GOTO_NUM_THREADS"] = "2"
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"
import datetime
import argparse
import glob
import copy

import asdf
from astropy.io import fits
from astropy.table import Table
from calcos import ccos

import numpy as np
import pandas as pd
import dask

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from query_cos_dark import files_by_mjd

class Superdark():
    def __init__(self, hv, segment, mjdstarts, mjdends, dayint=100, bin_x=1,
                 bin_y=1, bin_pha=1, phastart=1, phaend=31, pha_bins=None, 
                 gsagtab="/astro/sveash/cos_dark/5b919208l_gsag.fits",
                 bpixtab="/astro/sveash/cos_dark/5b91920bl_bpix.fits", 
                 region="inner", outfile=None,
                 outdir=".", xylimits=None):
        """
        To keep things consistent, start and stop range are defined the same
        as within python- start is inclusive but stop is EXCLUSIVE. E.g.
        using a PHA range of 1 to 31 means you will get PHA=1 up to PHA=30, but
        not PHA=31.
        """
        self.segment = segment
        self.hv = hv
        self.region = region
        self.mjdstarts = mjdstarts
        self.mjdends = mjdends
        self.ndays = [mjdends[i] - mjdstarts[i] for i in range(len(mjdstarts))]
        self.dayint = dayint
        self.bin_x = bin_x
        self.bin_y = bin_y
        
        self.phastart = phastart
        self.phaend = phaend
        if bin_pha == "all":
            bin_pha = phaend - phastart
        self.bin_pha = bin_pha
        self.pha_bins = pha_bins
        self.get_pha_bins()

        self.superdarks = []
        self.gsagtab = gsagtab
        self.bpixtab = bpixtab
        self.outfile = outfile
        self.outdir = outdir

        if xylimits == None:
            self.get_xy_limits()

        self.get_gsag_holes()
#        self.get_bpix_regions()


    @classmethod
    def from_asdf(cls, superdark):
        af = asdf.open(superdark)
        if "pha_bins" in af.keys():
            pha_bins = copy.deepcopy(af["pha_bins"])
        else:
            pha_bins = None
        self = cls(hv=af["hv"], segment=af["segment"], mjdstarts=af["mjdstarts"], 
                   mjdends=af["mjdends"], bin_x=af["bin_x"], bin_y=af["bin_y"],
                   bin_pha=af["bin_pha"], phastart=af["phastart"],
                   phaend=af["phaend"], pha_bins=pha_bins, 
                   outfile=superdark)
        self.outfile = superdark
        self.total_exptime = af["total_exptime"]
        self.total_files = af["total_files"]

        self.pha_images = {}
        for i in range(len(self.pha_bins)-1):
            key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
            self.pha_images[key] = copy.deepcopy(af[key])
            self.superdarks.append(copy.deepcopy(af[key]))
        af.close()
        return self


    def get_pha_bins(self):
        if self.pha_bins is None:
            self.pha_bins = np.arange(self.phastart, self.phaend+1, self.bin_pha)
        else:
            self.phastart = self.pha_bins[0]
            self.phaend = self.pha_bins[-1]
            grid = list(set(self.pha_bins[1:] - self.pha_bins[:-1]))
            if len(grid) == 1:
                self.bin_pha = grid[0]
            else:
                self.bin_pha = None


    def get_xy_limits(self):
        if self.segment == "FUVA":
            if self.region == "inner":
                xstart, xend, ystart, yend = 1260, 15119, 375, 660
        elif self.segment == "FUVB":
            if self.region == "inner":
                xstart, xend, ystart, yend = 1000, 14990, 405, 740
        else:
            raise Exception(f"Invalid segment specified: {segment}")
        while (xstart // self.bin_x) * self.bin_x < xstart:
            xstart += 1
        while (ystart // self.bin_y) * self.bin_y < ystart:
            ystart += 1
        self.xstart, self.ystart = xstart, ystart
        self.xend = (xend // self.bin_x) * self.bin_x
        self.yend = (yend // self.bin_y) * self.bin_y
        self.bin_xstart = self.xstart // self.bin_x
        self.bin_ystart = self.ystart // self.bin_y
        self.bin_xend = self.xend // self.bin_x
        self.bin_yend = self.yend // self.bin_y


    def create_superdark(self):

        @dask.delayed
        def _bin_dark(self, dark_df, phastart, phaend):
            pha_images = {} 
            darks = dark_df["fileloc"].values
            assert len(darks) != 0, "ERROR: no darks found at all!!"
                
            total_exptime = sum([fits.getval(x, "exptime", 1) for x in darks])
            gsag_df = self.gsag_holes.loc[(self.gsag_holes["DATE"] > self.mjdstarts[0]) & (self.gsag_holes["DATE"] < self.mjdends[0])]
#            bpix_df = self.bpix_regions
            print(f"   Binning corrtags, {phastart} <= PHA < {phaend}...")
            sum_image = self.bin_corrtags(darks, phastart=phastart, phaend=phaend)
            for j in range(len(gsag_df)):
                sum_image[gsag_df.iloc[j]["Y0"]:gsag_df.iloc[j]["Y1"], gsag_df.iloc[j]["X0"]:gsag_df.iloc[j]["X1"]] = 99999
#            for j in range(len(bpix_df)):
#                sum_image[bpix_df.iloc[j]["Y0"]:bpix_df.iloc[j]["Y1"], bpix_df.iloc[j]["X0"]:bpix_df.iloc[j]["X1"]] = 99999
            tmp = sum_image.reshape(1024 // self.bin_y, self.bin_y,
                                    16384 // self.bin_x, self.bin_x)
            binned = tmp.sum(axis=3).sum(axis=1)
            binned_inner = binned[self.bin_ystart:self.bin_yend, self.bin_xstart:self.bin_xend]
            zeros = np.where(binned_inner == 0)
            key = f"pha{phastart}-{phaend}"
            pha_images[key] = binned_inner
            return pha_images, total_exptime

        notfilled = True
        self.pha_images = {}
        runstart = datetime.datetime.now()
        print("\nStart time: {}".format(runstart))
        dark_dfs = []
        total_days = 0
        total_exptime = 0
        total_files = 0
        for k in range(len(self.mjdstarts)):
            start = self.mjdstarts[k]
            done = False
            while done is False:
                total_days += self.dayint
                if total_days > self.ndays[k]:
                    total_days = self.ndays[k]
                    done = True
                end = start + self.dayint
                if end > self.mjdends[k]:
                    end = self.mjdends[k]
                print(f"Using darks from MJD {start:,}-{end:,}; {total_days}/{self.ndays[k]} days")
                print("   Querying...")
                dark_df = files_by_mjd(start, end, segment=self.segment, hv=self.hv)
                dark_dfs.append(dark_df)
                print("   Query done")
                start = end

        dark_df = pd.concat(dark_dfs)
        delayed_bin = [_bin_dark(self, dark_df, self.pha_bins[i], self.pha_bins[i+1]) for i in range(len(self.pha_bins)-1)]
        out = dask.compute(*delayed_bin, scheduler='multiprocessing', num_workers=12)
        print("   Binning done")
        for grp in out:
            dct = grp[0]
            key = list(dct.keys())[0]
            self.pha_images[key] = dct[key]
        total_exptime += out[0][1]

        # Undo 99999 put in for gainsag holes
        for i in range(len(self.pha_bins)-1):
            key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
            inds = np.where(self.pha_images[key] >= 99999)
            print(f"Zeroing out {len(inds[0])} gain sagged pixels for {key}")
            self.pha_images[key][inds] = 0
            inds = np.where(self.pha_images[key] >= 99999)
            assert len(inds[0]) == 0, f"Not all gain sag holes were zeroed for {key}"
            self.superdarks.append(self.pha_images[key])

        self.total_exptime = total_exptime
        self.total_files = len(dark_df["fileloc"].values)

        # Write out ASDF file
        self.write_superdark(self.pha_images)
        runend = datetime.datetime.now()
        print("End time: {}".format(runend))
        print("Total time: {}".format(runend-runstart))

    
    def write_superdark(self, data_dict):
        data_dict["xstart"] = self.xstart
        data_dict["ystart"] = self.ystart
        data_dict["phastart"] = self.phastart
        data_dict["xend"] = self.xend
        data_dict["yend"] = self.yend
        data_dict["phaend"] = self.phaend
        data_dict["bin_x"] = self.bin_x
        data_dict["bin_y"] = self.bin_y
        data_dict["bin_pha"] = self.bin_pha
        data_dict["pha_bins"] = self.pha_bins
        data_dict["mjdstarts"] = self.mjdstarts
        data_dict["mjdends"] = self.mjdends
        data_dict["segment"] = self.segment
        data_dict["hv"] = self.hv
        data_dict["total_exptime"] = self.total_exptime
        data_dict["total_files"] = self.total_files
        af = asdf.AsdfFile(data_dict)
        today = datetime.datetime.now().strftime("%d%b%y-%H:%M:%S")
        if self.outfile is None:
            self.outfile = f"superdark_{self.segment}_{self.hv}_{today}.asdf"
        af.write_to(self.outfile)
        print(f"Wrote {self.outfile}")


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
    

    def get_bpix_regions(self):
        bpix = Table.read(self.bpixtab, format="fits", hdu=1)
        bpix_df = bpix.to_pandas()
        bpix_df = bpix_df.rename(columns={"LX": "X0", "LY": "Y0"})
        bpix_df["X1"] = bpix_df["X0"] + bpix_df["DX"]
        bpix_df["Y1"] = bpix_df["Y0"] + bpix_df["DY"]
        bpix_df["BIN_X0"] = bpix_df["X0"]//self.bin_x
        bpix_df["BIN_X1"] = bpix_df["X1"]//self.bin_x
        bpix_df["BIN_Y0"] = bpix_df["Y0"]//self.bin_y
        bpix_df["BIN_Y1"] = bpix_df["Y1"]//self.bin_y
        bpix_df = bpix_df.astype("int32")

        self.bpix_regions = bpix_df


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

    
    def bin_superdark(self, bin_x, bin_y, pha_bins=None, outfile=None, verbose=True):
        
        # Bin across PHA
        self.pha_images = {}
        if pha_bins is not None:
            for b in pha_bins:
                if b not in self.pha_bins:
                    raise IndexError(f"Previous PHA bins not compatible with new bins, {self.pha_bins} vs {pha_bins}")
            superdarks = []
            inds = np.nonzero(np.in1d(self.pha_bins, pha_bins))[0]
            for i in range(len(inds)-1):
                superdark = self.superdarks[i]
                for j in range(i+1, inds[i+1]):
                    superdark += self.superdarks[j]
                superdarks.append(superdark)
            self.pha_bins = np.array(pha_bins)
            self.get_pha_bins()
            self.superdarks = superdarks
       
        if outfile is None: 
            nowdt = datetime.datetime.now()
            now = nowdt.strftime("%d%b%Y")
            outfile = self.outfile.replace(".asdf", f"binned_{now}.asdf")
        self.outfile = outfile
        
        sh = np.shape(self.superdarks[0])
        xdim = sh[1]
        ydim = sh[0]
        b_x0 = 0
        b_y0 = 0
        b_x1 = (xdim // bin_x) * bin_x 
        b_y1 = (ydim // bin_y) * bin_y
        self.xend = b_x1 * self.bin_x + self.xstart
        self.yend = b_y1 * self.bin_y + self.ystart
        
        self.bin_x *= bin_x
        self.bin_y *= bin_y
        pdffile = os.path.join(self.outdir, self.outfile.replace("asdf", "pdf"))
        pdf = PdfPages(pdffile)
        # Bin in spatial directions
        for i,sd in enumerate(self.superdarks):
            phastart = self.pha_bins[i]
            phaend = self.pha_bins[i+1]

            binned = sd[b_y0:b_y1, b_x0:b_x1]
            binned_sh = np.shape(binned)
            binned_xdim = binned_sh[1]
            binned_ydim = binned_sh[0]
            tmp = binned.reshape(binned_ydim // bin_y, bin_y, binned_xdim // bin_x, bin_x)
            binned = tmp.sum(axis=3).sum(axis=1)
            self.superdarks[i] = binned
            key = f"pha{phastart}-{phaend}"
            self.pha_images[key] = binned
            rate = binned/self.total_exptime
            if verbose is True:
                print(f"For PHAs {phastart} through {phaend}")
                print(f"Binning by X={self.bin_x}, Y={self.bin_y}")
                print(f"\tTotal number of events: {np.sum(binned):,}")
                print(f"\tTotal exptime of superdark: {self.total_exptime:,}")
                print(f"\tMinimum number of events in a binned pixel: {np.min(binned)}")
                print(f"\tMaximum number of events in a binned pixel: {np.max(binned)}")
                print(f"\tMean number of events per binned pixel: {np.mean(binned):.1f}")
                print(f"\t  Standard deviation: {np.std(binned):.1f}")
                print(f"\tMedian number of events per binned pixel: {np.median(binned):.1f}")
                print(f"\tMean countrate per binned pixel: {np.mean(rate):.2e}")
                print(f"\t  Standard deviation: {np.std(rate):.2e}")
                print(f"\tMedian countrate per binned pixel: {np.median(rate):.2e}")
            fig, ax = plt.subplots(figsize=(20,5))
            #vmin = np.mean(rate) - 3*np.std(rate)
            vmin = np.median(rate) - np.median(rate)*0.5
            if vmin < 0:
                vmin = 0
            #vmax = np.mean(rate) + 3*np.std(rate)
            vmax = np.median(rate) + np.median(rate)*0.5
            im = ax.imshow(rate, aspect="auto",
                           origin="lower", cmap="inferno", vmin=vmin, vmax=vmax)
            fig.colorbar(im, label="Counts/s", format="%.2e")

            ax.set_title(f"{self.segment}; HV={self.hv}; MJD {self.mjdstarts}-{self.mjdends}; PHA {phastart}-{phaend}; X bin={self.bin_x} Y bin={self.bin_y}")
            plt.tight_layout()
            pdf.savefig(fig)
        pdf.close()
        print(f"Wrote {pdffile}")

        self.write_superdark(self.pha_images)

    def screen_counts(self, verbose=True, sigma=10):
        for i,sd in enumerate(self.superdarks):
            phastart = self.pha_bins[i]
            phaend = self.pha_bins[i+1]
            key = f"pha{phastart}-{phaend}"
            avg = np.mean(sd)
            std = np.std(sd)
            sigma_cutoff = (std*sigma)
            bad = np.where(sd > (avg+sigma_cutoff))
            if len(bad[0]) == 0:
                continue
            if verbose is True:
                print(f"{len(bad[0])} pixels have counts above {sigma}sigma, zeroing out now...")
            sd[bad] = 0.0
            self.superdarks[i] = sd
            self.pha_images[key] = sd

