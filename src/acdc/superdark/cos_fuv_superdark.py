import os
os.environ["GOTO_NUM_THREADS"] = "5"
os.environ["OMP_NUM_THREADS"] = "5"
os.environ["OPENBLAS_NUM_THREADS"] = "5"
os.environ["CRDS_SERVER_URL"] = "https://hst-crds.stsci.edu"
import crds
import datetime
import argparse
import glob
import copy
import asdf
from astropy.io import fits
from astropy.table import Table
from calcos import ccos
import numpy.ma as ma
import numpy as np
import pandas as pd
import dask
import matplotlib.pyplot as plt
dirname = os.path.dirname(__file__)
stylesheet = os.path.join(dirname, "../analysis", "niceplot.mplstyle")
plt.style.use(stylesheet)
from matplotlib import ticker
from matplotlib.backends.backend_pdf import PdfPages
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans

from acdc.database.query_cos_dark import files_by_mjd
from acdc.utils.utils import get_psa_wca, unbin_coords, bin_coords

class Superdark():
    def __init__(self, hv, segment, mjdstarts, mjdends, dayint=100, bin_x=1,
                 bin_y=1, bin_pha=1, phastart=1, phaend=31, pha_bins=None, 
                 gsagtab=None, bpixtab=None, region="inner", outfile=None,
                 outdir=".", xylimits=None, overwrite=False, sdqflags=8346):
                 #gsagtab="/grp/hst/cdbs/lref/71j1935gl_gsag.fits", 
                 #bpixtab="/grp/hst/cdbs/lref/69m20528l_bpix.fits", 
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
        self.overwrite = overwrite
        self.sdqflags = sdqflags

        self.phastart = phastart
        self.phaend = phaend
        if bin_pha == "all":
            bin_pha = phaend - phastart
        self.bin_pha = bin_pha
        self.pha_bins = pha_bins
        self.get_pha_bins()

        self.superdarks = []
        if gsagtab is None:
            gsagtab = self.get_reffile("gsagtab")
        self.gsagtab = gsagtab
        if bpixtab is None:
            bpixtab = self.get_reffile("bpixtab")
        self.bpixtab = bpixtab
        self.outfile = outfile
        self.outdir = outdir

        if xylimits == None:
            self.get_xy_limits()

        self.gsag_holes = get_gsag_holes(self.gsagtab, self.segment, self.hv)
        self.fixed_gsag = False 
        self.is_binned = False
#        self.get_bpix_regions()


    @classmethod
    def from_asdf(cls, superdark, overwrite=False):
        af = asdf.open(superdark)
        if "pha_bins" in af.keys():
            pha_bins = copy.deepcopy(af["pha_bins"])
        else:
            pha_bins = None
        self = cls(hv=af["hv"], segment=af["segment"], mjdstarts=af["mjdstarts"], 
                   mjdends=af["mjdends"], bin_x=af["bin_x"], bin_y=af["bin_y"],
                   bin_pha=af["bin_pha"], phastart=af["phastart"],
                   phaend=af["phaend"], pha_bins=pha_bins, xylimits="predetermined", 
                   outfile=superdark)
        for k in ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", 
                  "yend", "phastart", "phaend"]:
            setattr(self, k, af[k])

        self.outfile = superdark
        self.outfile_basename = os.path.basename(superdark)
        self.total_exptime = af["total_exptime"]
        self.total_files = af["total_files"]
        try:
            self.fixed_gsag = af["fixed_gsag"]
        except KeyError:
            self.fixed_gsag = False
        try:
            self.is_binned = af["is_binned"]
        except KeyError:
            if "binned" in os.path.basename(superdark):
                self.is_binned = True
            else:
                self.is_binned = False
        self.overwrite = overwrite

        self.pha_images = {}
        try:
            self.dq_image = copy.deepcopy(af["dq_image"])
        except KeyError:
            sh = af[f"pha{pha_bins[0]}-{pha_bins[1]}"].shape
            self.dq_image = np.zeros(sh[0]*sh[1]).reshape(sh).astype(int)
        for i in range(len(self.pha_bins)-1):
            key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
            self.pha_images[key] = copy.deepcopy(af[key])
            self.superdarks.append(copy.deepcopy(af[key]))

        try:
            self.creation_date = af["creation_date"]
        except KeyError:
            self.creation_date = "unknown"
        af.close()
        return self


    def get_reffile(self, filetype):
        now = datetime.datetime.now()
        current_pmap = crds.get_default_context()
        config = {"INSTRUME": "COS",
                  "CENWAVE": 1291,
                  "DETECTOR": "FUV",
                  "DATE-OBS": now.strftime("%Y-%m-%d"),
                  "TIME-OBS": "00:00:00"}
        crds_results = crds.getrecommendations(parameters=config, reftypes=[filetype],
                                               context=current_pmap, observatory="hst")
        reffile0 = crds_results[filetype]
        try:
            lref = os.environ["lref"]
        except KeyError as e:
            print(e.message)
            print("You must define the $lref environment variable- this is where all COS reference files are located")
        reffile = os.path.join(lref, reffile0)
        
        return reffile
    

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
            raise Exception(f"Invalid segment specified: {self.segment}")
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


    def create_superdark(self, dark_df=None):

        @dask.delayed
        def _bin_dark(self, dark_df, phastart, phaend):
            pha_images = {} 
            darks = dark_df["fileloc"].values
            assert len(darks) != 0, "ERROR: no darks found at all!!"
                
            total_exptime = sum([fits.getval(x, "exptime", 1) for x in darks])
            # Handle gainsag holes
            gsag_df = self.gsag_holes.loc[(self.gsag_holes["DATE"] > self.mjdstarts[0]) & 
                                          (self.gsag_holes["DATE"] < self.mjdends[0])] 
            print(f"   Binning corrtags, {phastart} <= PHA < {phaend}...")
            sum_image = self.bin_corrtags(darks, phastart=phastart, phaend=phaend)

            ylen,xlen = sum_image.shape
            dq_image = np.zeros((ylen, xlen)).astype(int)
            
            # Loop through each gainsag hole and mark the DQ image as DQ=1
            for j in range(len(gsag_df)):
                dq_image[gsag_df.iloc[j]["Y0"]:gsag_df.iloc[j]["Y1"], 
                    gsag_df.iloc[j]["X0"]:gsag_df.iloc[j]["X1"]] = 1

            ## Handle bad pixel regions
            #bpix_df = self.bpix_regions
            #bpix_df = bpix_df.loc[(bpix_df.SDQ == True) & (bpix_df.SEGMENT == self.segment)]
            ## Loop through each bad pixel region and set affected pixels to a false value of 99999
            #for j in range(len(bpix_df)):
            #    sum_image[bpix_df.iloc[j]["Y0"]:bpix_df.iloc[j]["Y1"], 
            #        bpix_df.iloc[j]["X0"]:bpix_df.iloc[j]["X1"]] = 99999

            ## Handle hotspots
            #spot_df = self.hotspots
            #spot_df = spot_df.loc[(spot_df.SEGMENT == self.segment) & (spot_df.START > self.mjdstarts[0]) & (spot_df.STOP < self.mjdends[0])]
            ## Loop through each hotspot and set affected pixels to a false value of 99999
            #for j in range(len(spot_df)):
            #    sum_image[spot_df.iloc[j]["Y0"]:spot_df.iloc[j]["Y1"], spot_df.iloc[j]["X0"]:spot_df.iloc[j]["X1"]] = 99999

            tmp = sum_image.reshape(1024 // self.bin_y, self.bin_y,
                                    16384 // self.bin_x, self.bin_x)
            binned = tmp.sum(axis=3).sum(axis=1)
            binned_inner = binned[self.bin_ystart:self.bin_yend, self.bin_xstart:self.bin_xend]
            tmp_dq = dq_image.reshape(1024 // self.bin_y, self.bin_y,
                                    16384 // self.bin_x, self.bin_x)
            binned_dq = tmp_dq.sum(axis=3).sum(axis=1)
            binned_inner_dq = binned_dq[self.bin_ystart:self.bin_yend, self.bin_xstart:self.bin_xend]
            key = f"pha{phastart}-{phaend}"
            pha_images[key] = binned_inner
            return pha_images, total_exptime, binned_inner_dq

        runstart = datetime.datetime.now()
        print("\nStart time: {}".format(runstart))
        total_exptime = 0
        total_files = 0
        self.pha_images = {}
        if dark_df is None:
            notfilled = True
            dark_dfs = []
            total_days = 0
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
        out = dask.compute(*delayed_bin, scheduler='multiprocessing', num_workers=15)
        print("   Binning done")
        for grp in out:
            dct = grp[0]
            dct_dq = grp[2]
            key = list(dct.keys())[0]
            self.pha_images[key] = dct[key]
        self.dq_image = out[0][2]
        total_exptime += out[0][1]

        ## Undo 99999 put in for gainsag holes
        #for i in range(len(self.pha_bins)-1):
        #    key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
        #    inds = np.where(self.pha_images[key] >= 99999)
        #    print(f"Zeroing out {len(inds[0])} gain sagged pixels for {key}")
        #    self.pha_images[key][inds] = 0
        #    inds = np.where(self.pha_images[key] >= 99999)
        #    assert len(inds[0]) == 0, f"Not all gain sag holes were zeroed for {key}"
        #    self.superdarks.append(self.pha_images[key])

        self.total_exptime = total_exptime
        self.total_files = len(dark_df["fileloc"].values)

        # Write out ASDF file
        self.write_superdark()
        runend = datetime.datetime.now()
        print("End time: {}".format(runend))
        print("Total time: {}".format(runend-runstart))

    
    def write_superdark(self, user_outfile=None, overwrite=False):
        data_dict = self.pha_images
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
        data_dict["fixed_gsag"] = self.fixed_gsag
        data_dict["is_binned"] = self.is_binned
        data_dict["dq_image"] = self.dq_image
        now = datetime.datetime.now()
        self.creation_date = now.strftime("%Y-%m-%d %H:%M")
        data_dict["creation_date"] = self.creation_date
        af = asdf.AsdfFile(data_dict)
        today = datetime.datetime.now().strftime("%d%b%y-%H:%M:%S")
        if self.outfile is None and user_outfile is None:
            self.outfile = f"superdark_{self.segment}_{self.hv}_{today}.asdf"
        elif self.outfile is None and user_outfile is not None:
            self.outfile = user_outfile
        if os.path.exists(self.outfile):
            assert sum([overwrite, self.overwrite]) > 0, f"Output superdark {self.outfile} already exists and overwrite is False"

        af.write_to(self.outfile)
        print(f"Wrote {self.outfile}")


    def fix_gsag(self, method="interpolate", kernel=None, row_threshold=.3,
                 lp1_interpolate=None):
        superdarks = []

        if self.segment == "FUVA":
            _, lp1_lims = bin_coords([0, 0], [459, 528], self.bin_x,
                self.bin_y, self.xstart, self.ystart)
        else:
            _, lp1_lims = bin_coords([0, 0], [521, 587], self.bin_x,
                self.bin_y, self.xstart, self.ystart)

        print(f"Looking for rows where >= {int(row_threshold*100)}% of the row is gain-sagged...")
        for i in range(len(self.pha_bins)-1):
            key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
            ylen,xlen = self.pha_images[key].shape
            for j in range(ylen):
                #if lp1_lims[0] <= j <= lp1_lims[1] and lp1_interpolate is True:
                #    print(f"    lp1_interpolate is True, interpolating over row {j}")
                #    self.dq_image[j,:] = 1
                #    continue
                if lp1_lims[0] <= j <= lp1_lims[1] and lp1_interpolate is False:
                    continue
                sagged = np.where(self.dq_image[j,:] == 1)
                frac_sagged = len(sagged[0])/xlen
                if frac_sagged >= row_threshold:
                    print(f"    {frac_sagged*100:.1f}% of row {j} is gain-sagged, entire row will be corrected")
                    self.dq_image[j,:] = 1

        inds = np.where(self.dq_image == 1)
        # Zero out gain sagged regions
        if method == "zero":
            print(f"Zeroing out gain sagged pixels...")
            for i in range(len(self.pha_bins)-1):
                key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
                self.pha_images[key][inds] = 0
                #inds = np.where(self.pha_images[key] >= 99999)
                #assert len(inds[0]) == 0, f"Not all gain sag holes were zeroed for {key}"
                superdarks.append(self.pha_images[key])
        elif method == "boost":
            print("Boosting gain sagged pixels...")
            for i in range(len(self.pha_bins)-1):
                key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
                indsbelow = (inds[0]-1, inds[1])
                indsabove = (inds[0]+1, inds[1])
                avg_around = (self.pha_images[key][indsbelow]+self.pha_images[key][indsabove])/2
                self.pha_images[key][inds] = avg_around
                superdarks.append(self.pha_images[key])
        elif method == "interpolate":
            print("Interpolating over gain sagged pixels...")
            if kernel is None:
                kernel = Gaussian2DKernel(x_stddev=.3, y_stddev=.3, x_size=11, y_size=5)
            for i in range(len(self.pha_bins)-1):
                key = f"pha{self.pha_bins[i]}-{self.pha_bins[i+1]}"
                self.pha_images[key][inds] = np.nan
                interp_sd = interpolate_replace_nans(self.pha_images[key], kernel)
                self.pha_images[key] = interp_sd
                superdarks.append(self.pha_images[key])
        self.superdarks = superdarks
        self.fixed_gsag = True 


    def get_hotspots(self, extnum):
        spot = Table.read(self.spottab, format="fits", hdu=extnum)
        spot_df = spot.to_pandas()
        spot_df = spot_df.rename(columns={"LX": "X0", "LY": "Y0"})
        spot_df["X1"] = spot_df["X0"] + spot_df["DX"]
        spot_df["Y1"] = spot_df["Y0"] + spot_df["DY"]
        spot_df = spot_df.astype("int32")
        str_df = spot_df.select_dtypes([np.object])
        str_df = str_df.stack().str.decode('utf-8').unstack()
        for col in str_df:
            spot_df[col] = str_df[col]

        self.hotspots = spot_df
    

    def get_bpix_regions(self):
        bpix = Table.read(self.bpixtab, format="fits", hdu=1)
        bpix_df = bpix.to_pandas()
        bpix_df = bpix_df.rename(columns={"LX": "X0", "LY": "Y0"})
        bpix_df["X1"] = bpix_df["X0"] + bpix_df["DX"]
        bpix_df["Y1"] = bpix_df["Y0"] + bpix_df["DY"]
        # This is to fix the astropy bytes string columns 
        str_df = bpix_df.select_dtypes([np.object])
        str_df = str_df.stack().str.decode('utf-8').unstack()
        for col in str_df:
            bpix_df[col] = str_df[col]
        bpix_df["SDQ"] =  np.where(bpix_df["DQ"]&self.sdqflags == bpix_df["DQ"], True, False)

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
                print(f"No counts for PHA>={phastart} and PHA<{phaend}")
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


    def bin_superdark(self, bin_x, bin_y, pha_bins=None, outfile=None, verbose=True, writefile=True, makeplots=True):
        
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
        
        infile = self.outfile
        if outfile is None and writefile is True:
            nowdt = datetime.datetime.now()
            now = nowdt.strftime("%d%b%Y")
            outfile = self.outfile.replace(".asdf", f"_binned_{now}.asdf")
        self.outfile = outfile
        if writefile is True and self.overwrite is False and os.path.exists(outfile):
            print(f"WARNING: Output superdark {outfile} already exists and overwrite is False, skipping...")
            return

        sh = np.shape(self.superdarks[0])
        xdim = sh[1]
        ydim = sh[0]
        b_x0 = 0
        b_y0 = 0
        b_x1 = (xdim // bin_x) * bin_x 
        b_y1 = (ydim // bin_y) * bin_y
        self.xend = b_x1 * self.bin_x + self.xstart
        self.yend = b_y1 * self.bin_y + self.ystart
# TODO
        # Should this really be *=, not just =?        
        self.bin_x *= bin_x
        self.bin_y *= bin_y
        # Bin DQ image
        binned_dq = self.dq_image[b_y0:b_y1, b_x0:b_x1]
        binned_sh_dq = np.shape(binned_dq)
        binned_xdim_dq = binned_sh_dq[1]
        binned_ydim_dq = binned_sh_dq[0]
        tmp_dq = binned_dq.reshape(binned_ydim_dq // bin_y, bin_y, binned_xdim_dq // bin_x, bin_x)
        tmp_dq2 = np.bitwise_or.reduce(tmp_dq, axis=3)
        binned_dq = np.bitwise_or.reduce(tmp_dq2, axis=1)
        self.dq_image = binned_dq

        # Bin superdark in spatial directions
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
            zeroinds = np.where(binned == 0)
            nzero = len(zeroinds[0])
            if verbose is True:
                print("")
                print(f"Binning superdark {infile}")
                print(f"Binning PHAs {phastart} through {phaend}")
                print(f"Binning by X={self.bin_x}, Y={self.bin_y}")
                print(f"\tTotal number of events: {np.sum(binned):,}")
                print(f"\tTotal exptime [s] of superdark: {self.total_exptime:,.1f}")
                print(f"\tMinimum number of events in a binned pixel: {np.min(binned)}")
                print(f"\tMaximum number of events in a binned pixel: {np.max(binned)}")
                print(f"\tMean number of events per binned pixel: {np.mean(binned):.1f}")
                print(f"\t  Standard deviation: {np.std(binned):.1f}")
                print(f"\tMedian number of events per binned pixel: {np.median(binned):.1f}")
                print(f"\tMean countrate per binned pixel: {np.mean(rate):.2e}")
                print(f"\t  Standard deviation: {np.std(rate):.2e}")
                print(f"\tMedian countrate per binned pixel: {np.median(rate):.2e}")
                print(f"\tNumber of binned pixels with zero events: {nzero:,}")

        self.is_binned = True 

        if writefile is True:
            self.write_superdark()
        if makeplots is True: 
            self.plot_superdarks()

    def plot_superdarks(self, pdffile=None, vmin=None, vmax=None, savepng=False, pngroot=None):
        if savepng is False:
            if pdffile is None:
                pdffile = os.path.join(self.outdir, self.outfile.replace("asdf", "pdf"))
            pdf = PdfPages(pdffile)
        if savepng is True and pngroot is None:
            pngroot = os.path.join(self.outdir, self.outfile.replace(".asdf", ""))
        for i,sd in enumerate(self.superdarks):
            phastart = self.pha_bins[i]
            phaend = self.pha_bins[i+1]
            rate = sd/self.total_exptime
            fig, ax = plt.subplots(figsize=(20,5))
            #vmin = np.mean(rate) - 3*np.std(rate)
            if vmin is None:
                vmin = np.median(rate) - np.median(rate)*0.5
                if vmin < 0:
                    vmin = 0
            #vmax = np.mean(rate) + 3*np.std(rate)
            if vmax is None:
                vmax = np.median(rate) + np.median(rate)*1.
            im = ax.imshow(rate, aspect="auto", interpolation="nearest", 
                           origin="lower", cmap="inferno", vmin=vmin, vmax=vmax)
            cbar = fig.colorbar(im, label="Counts/s", format="%.1e", pad=0.01)
            cbarticks = cbar.get_ticks()
            # If there are more than 5 tickmarks, reduce it.
            if len(cbarticks) > 7: # 1st and last ticks are not shown
                cbar.ax.locator_params(nbins=5)
            ax.set_title(f"{self.segment}; HV={self.hv}; MJD {self.mjdstarts}-{self.mjdends}; PHA {phastart}-{phaend}; X bin={self.bin_x}, Y bin={self.bin_y}")
            ax.set_xlabel("X (binned)")
            ax.set_ylabel("Y (binned)")
            xticks = ax.xaxis.get_major_ticks()
            xticks[0].set_visible(False)
            xticks[1].set_visible(False)
            if savepng is False:
                pdf.savefig(fig, bbox_inches="tight")
            else:
                pngfile = pngroot+f"_pha{phastart}-{phaend}.png"
                fig.savefig(pngfile, bbox_inches="tight", dpi=200)
                print(f"Wrote {pngfile}")
        if savepng is False:
            pdf.close()
            print(f"Wrote {pdffile}")
        plt.close('all')
        plt.close(fig)


    # Should investigate if this can be replaced with typical sigma clipping
    def screen_counts(self, verbose=True, sigma=10, mask=False, interpolate=False, interp_kernel=None,  method=np.median, exclude_zeros=True):
        for i,sd in enumerate(self.superdarks):
            phastart = self.pha_bins[i]
            phaend = self.pha_bins[i+1]
            key = f"pha{phastart}-{phaend}"
            if exclude_zeros is True:
                nonzeroinds = np.where(sd > 0)
                sd_stats = sd[nonzeroinds]
            else:
                sd_stats = sd
            mid = method(sd_stats)
            std = np.std(sd_stats)
            sigma_cutoff = (std*sigma)
            bad = np.where(sd > (mid+sigma_cutoff))
            if len(bad[0]) == 0:
                continue
            if verbose is True:
                print(f"{len(bad[0])} pixels have counts above {sigma}sigma for {key}")

            if interpolate is True:
                if interp_kernel is None:
                    kernel = Gaussian2DKernel(x_stddev=.3, y_stddev=.3, x_size=11, y_size=5)
                # From https://docs.astropy.org/en/stable/convolution/index.html
                sd[bad] = np.nan
                interp_sd = interpolate_replace_nans(sd, kernel)
                self.superdarks[i] = interp_sd
                self.pha_images[key] = interp_sd
            elif mask is True:
                sd_masked = ma.masked_greater_equal(sd, mid+sigma_cutoff)
                self.superdarks[i] = sd_masked
                self.pha_images[key] = sd_masked
            else:
                sd[bad] = 0.0
                self.superdarks[i] = sd
                self.pha_images[key] = sd


    def write_fits(self, fitsfile=None):
        if fitsfile is None:
            nowdt = datetime.datetime.now()
            now = nowdt.strftime("%d%b%Y")
            fitsfile = self.outfile.replace(".asdf", ".fits")
        self.fitsfile = fitsfile
        if self.overwrite is False and os.path.exists(fitsfile):
            print(f"{fitsfile} already exists and overwrite is False")
            return
        hdr0 = fits.Header()
        hdr0["suprdark"] = os.path.basename(self.outfile)
        hdr0["hv"] = self.hv
        hdr0["segment"] = self.segment 
        primary = fits.PrimaryHDU(header=hdr0)
        hduims = []
        for i,sd in enumerate(self.superdarks):
            hdr = fits.Header()
            hdr["EXTNAME"] = "SCI"
            hdr["phastart"] = self.pha_bins[i]
            hdr["phaend"] = self.pha_bins[i+1]
            hdr["xstart"] = self.xstart
            hdr["ystart"] = self.ystart
            hdr["xend"] = self.xend
            hdr["yend"] = self.yend
            hdr["xbinsize"] = self.bin_x
            hdr["ybinsize"] = self.bin_y
            hdu = fits.ImageHDU(sd, header=hdr)
            hduims.append(hdu)
        hdulist = fits.HDUList([primary] + hduims)
        hdulist.writeto(fitsfile, overwrite=self.overwrite)
        print(f"Wrote {fitsfile}") 


    def write_cube_fits(self, fitsfile=None):
        if fitsfile is None:
            nowdt = datetime.datetime.now()
            now = nowdt.strftime("%d%b%Y")
            fitsfile = self.outfile.replace(".asdf", "_cube.fits")
        self.fitsfile = fitsfile
        if self.overwrite is False and os.path.exists(fitsfile):
            print(f"{fitsfile} already exists and overwrite is False")
            return
        hdr0 = fits.Header()
        hdr0["suprdark"] = os.path.basename(self.outfile)
        hdr0["hv"] = self.hv
        hdr0["segment"] = self.segment 
        primary = fits.PrimaryHDU(header=hdr0)
        cube = np.stack(self.superdarks, axis=0)
        hdr1 = fits.Header()
        hdr1["EXTNAME"] = "SCI"
        hdr1["phastart"] = self.pha_bins[0]
        hdr1["phaend"] = self.pha_bins[-1]
        hdr1["xstart"] = self.xstart
        hdr1["ystart"] = self.ystart
        hdr1["xend"] = self.xend
        hdr1["yend"] = self.yend
        hdr1["xbinsize"] = self.bin_x
        hdr1["ybinsize"] = self.bin_y
        hdu1 = fits.ImageHDU(cube, header=hdr1)
        hdulist = fits.HDUList([primary, hdu1])
        hdulist.writeto(fitsfile, overwrite=self.overwrite)
        print(f"Wrote {fitsfile}") 


    def get_binning_pars(self):
        keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
                "phastart", "phaend"]
        binning = {}
        for k in keys:
            binning[k] = getattr(self, k)

        return binning


def get_gsag_holes(gsagtab, segment, hv):
    extfound = False
    extnum = 1
    hvkey = f"HVLEVEL{segment[-1]}"
    while extfound is False:
        try:
            gsag_hv = fits.getval(gsagtab, hvkey, extnum)
            gsag_seg = fits.getval(gsagtab, "segment", extnum)
        except KeyError:
            extnum += 1
            continue
        if gsag_hv == hv and gsag_seg == segment:
            extfound = True
        else:
            extnum += 1
    gsag = Table.read(gsagtab, format="fits", hdu=extnum)
    gsag_df = gsag.to_pandas()
    gsag_df = gsag_df.rename(columns={"LX": "X0", "LY": "Y0"})
    gsag_df["X1"] = gsag_df["X0"] + gsag_df["DX"]
    gsag_df["Y1"] = gsag_df["Y0"] + gsag_df["DY"]
    gsag_df = gsag_df.astype("int32")

    return gsag_df
