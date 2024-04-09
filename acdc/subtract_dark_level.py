"""
Perform a custom dark subtraction on COS corrtag files.
"""

import os
import argparse
import copy
import asdf
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
dirname = os.path.dirname(__file__)
stylesheet = os.path.join(dirname, "analysis", "niceplot.mplstyle")
plt.style.use(stylesheet)
import glob

from .predict_dark_level import bin_science, get_binning_pars

RESEL = [6, 10]
PHA_INCLUSIVE = [2, 23]
PHA_INCL_EXCL = [2, 24]

def subtract_dark(corrtags, datadir, fact=1, outdir=".", overwrite=False):
    """Perform the custom dark correction.

    Args:
        corrtags (list or array-like): Input, unmodified science corrtags.
        datadir (str): Directory that houses the modeled dark rates for each
            input corrtag.
        fact (int): Changes resolution, 1=highest resolution, 8=lowest resolution. 
            To keep at same input resolution, use 1.
        outdir (str): Directory where custom corrtags will be written.
        """

    for item in corrtags:
        try:
            acdc_done = fits.getval(item, "ACDCCORR")
            if acdc_done == "COMPLETE":
                print(f"WARNING:\nCustom dark correction has already been applied to {item}, skipping...")
                continue
        except KeyError:
            pass
        
        outfile = os.path.join(outdir, f"corrected_{os.path.basename(item)}")
        if os.path.exists(outfile) and overwrite is False:
             print(f"WARNING: Corrected corrtag {outfile} already exists and overwrite is False, skipping...")
             continue
        if not os.path.exists(outdir):
            os.makedirs(outdir)
       
        hdr0 = fits.getheader(item, 0) 
        rootname = hdr0["rootname"]
        segment = hdr0["segment"]
        cenwave = hdr0["cenwave"]
        lp = hdr0["life_adj"]
        pred_noise_file = os.path.join(datadir, f"{rootname}_{segment}_noise_complete.asdf")
        pred_noise_af = asdf.open(pred_noise_file)
        pha_str = f"pha{PHA_INCL_EXCL[0]}-{PHA_INCL_EXCL[1]}"
        pred_noise = pred_noise_af[pha_str]
        b = get_binning_pars(pred_noise_af)
        binned, nevents = bin_science(item, b, segment, cenwave, lp, fact, exclude_airglow=False)
        
        hdulist = fits.open(item)
        data = hdulist[1].data
        
        xstart = b["xstart"]
        xend = b["xend"]
        ystart = b["ystart"]
        yend = b["yend"]
        phastart = b["phastart"]
        phaend = b["phaend"]
        bin_x = b["bin_x"]
        bin_y = b["bin_y"]

        #wwf = open("tst.txt", "w")
        logic = np.zeros((binned.shape[0], binned.shape[1]))
       
        all_delta_eps_ind = []
        for i in range(len(data["xcorr"])):
            if data["xcorr"][i] > xstart and data["xcorr"][i] < xend:
                if data["ycorr"][i] > ystart and data["ycorr"][i] < yend:
                    if data["pha"][i] > phastart and data["pha"][i] < phaend:
                        x = int((data["xcorr"][i] - xstart) // bin_x * fact)
                        y = int((data["ycorr"][i] - ystart) // bin_y)
                        delta_eps = pred_noise[y, int(x/fact)] / float(fact)
                        logic[y,x] = 1
                        delta_eps_ind = delta_eps / float(nevents[y,x])
                        if np.isnan(delta_eps_ind):
                            delta_eps_ind = 0
        #                outstr = f"{data['epsilon'][i]}\t{delta_eps}\t{delta_eps_ind}\n"
                        data["epsilon"][i] -= delta_eps_ind
        #                wwf.write(outstr)
        #wwf.close()

        #plt.imshow(logic, aspect="auto", origin="lower")
        #plt.savefig("logic.png")
        #plt.clf()
        
#        hdulist = fits.open(item, mode="update")
#        data = hdulist[1].data
#        phainds = np.where((data["pha"] >= b["phastart"]) &
#                           (data["pha"] <= b["phaend"]))
#        phadata = data[phainds]
#        innerinds = np.where((phadata["xcorr"] > b["xstart"]) &
#                             (phadata["xcorr"] < b["xend"]) &
#                             (phadata["ycorr"] > b["ystart"]) &
#                             (phadata["ycorr"] < b["yend"]))
#        filtered = phadata[innerinds]
#
#        wwf = open("tst.txt", "w")
#        logic = np.zeros((binned.shape[0], binned.shape[1]))
##        xs = (filtered["xcorr"] - b["xstart"]) // b["bin_x"] * fact
##        xs = (filtered["ycorr"] - b["ystart"]) // b["bin_y"]
##        delta_eps = pred_noise[
#        print("before loop 1")
#        all_delta_eps_ind = []
#        for i in range(len(filtered["xcorr"])):
#            print(i)
#            x = int((filtered["xcorr"][i] - b["xstart"]) // b["bin_x"] * fact)
#            y = int((filtered["ycorr"][i] - b["ystart"]) // b["bin_y"])
#            delta_eps = pred_noise[y, int(x/fact)] / float(fact)
#            logic[y,x] = 1
#            delta_eps_ind = delta_eps / float(nevents[y,x])
##            all_delta_eps_ind.append(delta_eps_ind)
#            outstr = f"{data[phainds][innerinds][i]}\t{delta_eps}\t{delta_eps_ind}\n"
#        #    data[phainds][innerinds]["epsilon"][i] -= delta_eps_ind
#            filtered["epsilon"][i] -= delta_eps_ind
#            wwf.write(outstr)
#        
#        data[phainds][innerinds]["epsilon"] = all_delta_eps_ind
##        data[phainds][innerinds]["epsilon"] -= all_delta_eps_ind
#        print("after loop 1")
#        wwf.close()
#
#        plt.imshow(logic, aspect="auto", origin="lower")
#        plt.savefig("logic.png")
#        plt.clf()

        
        inds = np.where(logic < 1)
        if len(logic[inds]) != 0:
            if len(logic[inds]) >  len(data):
                #num = int(np.ceil(len(logic[inds]) / len(data)))
                #tmp = [data]*num
                #longerdata = np.concatenate(tmp)
                longerhdu = fits.BinTableHDU.from_columns(hdulist[1].columns, nrows=len(logic[inds]))
                for colname in hdulist[1].columns.names:
                    longercol = np.resize(data[colname], len(logic[inds]))
                    longerhdu.data[colname][:] = longercol
                new_events = longerhdu.data
            else:
                new_events = copy.deepcopy(data[len(data) - len(inds[0]):])
            for i in range(len(inds[0])):
                y = b["bin_y"] * inds[0][i] + b["ystart"] 
                x = b["bin_x"] * inds[1][i] + b["xstart"]
                new_events["rawx"][i] = x 
                new_events["rawy"][i] = y
                new_events["xcorr"][i] = x 
                new_events["ycorr"][i] = y
                new_events["xfull"][i] = x 
                new_events["yfull"][i] = y
                new_events["pha"][i] = 10
                new_events["dq"][i] = 0
                new_events["epsilon"][i] = -pred_noise[inds[0][i], int(inds[1][i]/fact)] / fact
            data = np.concatenate([data, new_events])
        nans = np.isnan(data["EPSILON"])
        data["EPSILON"][nans] = 0
        hdulist[1].data = data
        hdulist.writeto(outfile, overwrite=overwrite)
        with fits.open(outfile, mode="update") as hdulist:
            hdr0 = hdulist[0].header
            hdr0.set("ACDCCORR", "COMPLETE")
        print(f"Wrote {outfile}")

        pred_noise_af.close() 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datadir",
                        help="Path to science corrtags")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Output directory")
    args = parser.parse_args()
    corrtags = glob.glob(os.path.join(args.datadir, "*corrtag*fits"))

    subtract_dark(corrtags, args.datadir, outdir=args.outdir)

