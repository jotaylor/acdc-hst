"""
Given COS science datasets, compare the background regions in the science
exposures to a low/quiescent and high/active superdark.

Command line-arguments:
-d or --datatdir
    The path to science corrtags
--hv 
    If datadir contains corrtags of multiple HVs, you can use --hv to use
    just one specific HV value.
--segment
    If datadir contains corrtags of multiple segments, you can use --segment
    to use just one specific segment.
--lo
    The name of the low activity superdark.
--hi
    The name of the high activity superdark.
-o or --outdir
    The name of the output directory where products will be written.
--binned
    Using this argument indicates that the supplied superdarks are already
    binned and should not be binned any further.
"""

import copy
import argparse
import os
import glob
import asdf
from scipy.optimize import minimize
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

from acdc.superdark.cos_fuv_superdark import Superdark
from acdc.database.calculate_dark import get_aperture_region
from acdc.analysis.compare_backgrounds import smooth_array
from acdc.utils.utils import bin_coords

RESEL = [6, 10]
PHA_INCLUSIVE = [2, 23]
PHA_INCL_EXCL = [2, 24]
# Colorblind-safe palette below
COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6cee3", 
          "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", 
          "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#a6cee3"]

def fun_opt(coeffs, darks, binned_sci, excluded_rows):
	combined_dark = linear_combination(darks, coeffs)
	cval = c_stat(combined_dark, binned_sci, excluded_rows)
	
	return cval


def linear_combination(darks, coeffs):
    combined = darks[0] * coeffs[0]
    for i in range(1, len(darks)):
        combined += darks[i] * coeffs[i]
    return combined 


def check_superdarks(af1, af2):
    """
    Ensure that both active and quiescent superdarks were binned in identical
    ways.
    """
    bad = False
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
            "phastart", "phaend"]
    for k in keys:
        if af1[k] != af2[k]:
            print(f"WARNING!!!! Key {k} does not match for both superdarks")
            bad = True
    assert bad == False, "Cannot continue until discrepancies are resolved"


def get_binning_pars(af):
    """
    For a given superdark, return the binning information in the spatial 
    directions and PHA.
    """
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
            "phastart", "phaend"]
    binning = {}
    for k in keys:
        binning[k] = af[k]

    return binning

def bin_science(corrtag, b, fact=1):
    """
    Given a corrtag with lists of events as a function of X, Y, and PHA,
    bin the data into an image using the same bin sizes as the superdark.
    """
    data = fits.getdata(corrtag)
    phainds = np.where((data["pha"] >= b["phastart"]) & 
                       (data["pha"] < b["phaend"]))
    phadata = data[phainds]
    innerinds = np.where((phadata["xcorr"] > b["xstart"]) &
                         (phadata["xcorr"] < b["xend"]) & 
                         (phadata["ycorr"] > b["ystart"]) &
                         (phadata["ycorr"] < b["yend"]))
    filtered = phadata[innerinds]

    ydim = int((b["yend"] - b["ystart"]) / b["bin_y"])
    xdim = int((b["xend"] - b["xstart"]) / b["bin_x"])
    binned = np.zeros((ydim, xdim*fact))
    nevents = np.zeros((ydim, xdim*fact))
    
    for i in range(len(filtered["xcorr"])):
        x = int((filtered["xcorr"][i] - b["xstart"]) // b["bin_x"] * fact)
        y = int((filtered["ycorr"][i] - b["ystart"]) // b["bin_y"])
        binned[y,x] += filtered["epsilon"][i]
        nevents[y,x] += 1
    return binned, nevents


def get_excluded_rows(segment, cenwave, lp, binning):
    """
    Determine the rows that correspond to the PSA and WCA apertures for a 
    given segment and lifetime position. Return the row indices in the
    binned superdark coordinate system.
    """
    
    excluded_rows = np.array(())
    apertures = {"PSA": None, "WCA": None}
    for aperture in apertures:
        aperture_regions = get_aperture_region(cenwave=cenwave, aperture=aperture)
        box = aperture_regions[segment][f"lp{lp}_{aperture.lower()}_{cenwave}"]
        xmin0, xmax0, ymin0, ymax0 = box
        xsnew, ysnew = bin_coords(np.array([xmin0, xmax0]), np.array([ymin0, ymax0]), binning["bin_x"], binning["bin_y"], binning["xstart"], binning["ystart"], make_int=True)
        ap_xmin, ap_xmax = xsnew
        ap_ymin, ap_ymax = ysnew
        rows = np.arange(ap_ymin, ap_ymax+1)
        apertures[aperture] = rows
        excluded_rows = np.concatenate((excluded_rows, rows))
    excluded_rows = excluded_rows.astype(int)
    return excluded_rows, apertures

def c_stat(combined_dark, binned_sci, excluded_rows):
    csum = 0.
    for y in range(combined_dark.shape[0]):
        for x in range(combined_dark.shape[1]):
            if binned_sci[y,x] > 0 and y not in excluded_rows and combined_dark[y,x] > 0:
                csum = csum + 2. * (combined_dark[y,x] - binned_sci[y,x] + binned_sci[y,x] * (np.log(binned_sci[y,x]) - np.log(combined_dark[y,x])))
            elif y not in excluded_rows:
                csum += np.abs(combined_dark[y,x]) * 2.
            elif combined_dark[y,x] <= 0:
                csum += np.inf

    return csum


def predict_dark(corrtags, lo_darkname, hi_darkname, segment=None, hv=None, outdir=".", binned=False):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if binned is False:
        Lo = Superdark.from_asdf(lo_darkname)
        Hi = Superdark.from_asdf(hi_darkname)
        lo_binnedname = lo_darkname.replace(".asdf", "_phabinned.asdf")
        hi_binnedname = hi_darkname.replace(".asdf", "_phabinned.asdf")
        Lo.screen_counts(verbose=False)
        Hi.screen_counts(verbose=False)
        Lo.bin_superdark(RESEL[0]*2, RESEL[1]*2, pha_bins=PHA_INCL_EXCL, outfile=lo_binnedname)
        Hi.bin_superdark(RESEL[0]*2, RESEL[1]*2, pha_bins=PHA_INCL_EXCL, outfile=hi_binnedname)
    else:
        lo_binnedname = lo_darkname
        hi_binnedname = hi_darkname
    lo_af = asdf.open(lo_binnedname)
    hi_af = asdf.open(hi_binnedname)
    check_superdarks(lo_af, hi_af)
    binning = get_binning_pars(lo_af)
    pha_str = f"pha{PHA_INCL_EXCL[0]}-{PHA_INCL_EXCL[1]}"
    lo_dark = lo_af[pha_str]
    hi_dark = hi_af[pha_str]
#    lo_dark = lo_dark[:, :288]
#    hi_dark = hi_dark[:, :288]
    dark_hv = lo_af["hv"]
    dark_segment = lo_af["segment"]
    fig,ax = plt.subplots(1, 1, figsize=(20,8))
    im = ax.imshow(lo_dark, aspect="auto", origin="lower")
    fig.colorbar(im, label="Counts")
    ax.set_title(f"{dark_segment} {dark_hv} quiescent superdark", size=25)
    figname = os.path.join(outdir, f"{dark_segment}_{dark_hv}_lo_superdark.png")
    fig.savefig(figname, bbox_inches="tight")
    print(f"Saved quiscent superdark figure: {figname}")
    fig,ax = plt.subplots(1, 1, figsize=(20,8))
    im = ax.imshow(hi_dark, aspect="auto", origin="lower")
    fig.colorbar(im, label="Counts")
    ax.set_title(f"{dark_segment} {dark_hv} active superdark", size=25)
    figname = os.path.join(outdir, f"{dark_segment}_{dark_hv}_hi_superdark.png")
    fig.savefig(figname, bbox_inches="tight")
    print(f"Saved active superdark figure:{figname}")
    for item in corrtags:
        if segment is not None:
            file_segment = fits.getval(item, "segment")
            if hv is not None:
                file_hv = fits.getval(item, f"HVLEVEL{file_segment[-1]}", 1)
                if file_segment != segment.upper() and str(file_hv) != str(hv):
                    print(f"File does not match required settings: {item}")
                    continue
        else:
            segment = fits.getval(item, "segment") 

        cenwave = fits.getval(item, "cenwave")
        fig, ax = plt.subplots(1, 1, figsize=(20, 8))
        lp = fits.getval(item, "life_adj")
        excluded_rows, apertures = get_excluded_rows(segment, cenwave, lp, binning)
        rootname0 = fits.getval(item, "rootname")
        rootname = rootname0.lower()
        binned_sci, nevents = bin_science(item, binning)
        vmin = np.mean(binned_sci) - 1*np.std(binned_sci)
        if vmin < 0:
            vmin = 0
        vmax = np.mean(binned_sci) + 1*np.std(binned_sci)
        sh = np.shape(binned_sci)
        im = ax.imshow(binned_sci, aspect="auto", origin="lower", 
                vmin=vmin, vmax=vmax, extent=[0, sh[1], 0, sh[0]])
        ax.axhline(apertures["PSA"][0], lw=2, color=COLORS[3], label="PSA")
        ax.axhline(apertures["PSA"][-1]+1, lw=2, color=COLORS[3])
        ax.axhline(apertures["WCA"][0], lw=2, ls="dashed", color=COLORS[3], label="WCA")
        ax.axhline(apertures["WCA"][-1]+1, lw=2, ls="dashed", color=COLORS[3])
        ax.legend(loc="upper right")
        fig.colorbar(im, label="Counts")
        figname = os.path.join(outdir, f"{rootname}_{segment}_binned_sci.png")
        ax.set_title(f"{rootname} Binned Science Image", size=25)
        fig.savefig(figname, bbox_inches="tight")
        print(f"Saved binned science image: {figname}")
        plt.clf()
        sci_exp = fits.getval(item, "exptime", 1)
        lo_exptime = lo_af["total_exptime"]
        hi_exptime = hi_af["total_exptime"]
        # Initial guess is 0.5 contribution for each superdark
        lo_coeff = 0.5 / (lo_exptime / sci_exp)
        hi_coeff = 0.5 / (hi_exptime / sci_exp)
        combined_dark = linear_combination([lo_dark, hi_dark], [lo_coeff, hi_coeff])
# this is science extraction regions (use xtractab) for sci exp. LP
# also exclude wavecals
        nrows = np.shape(binned_sci)[0]
        all_rows = np.arange(nrows)
        bkg_rows = list(set(all_rows) - set(excluded_rows))
        plt.figure(figsize=(20, 8))
        for i in range(binned_sci.shape[0]):
            if i not in excluded_rows:
                plt.plot(binned_sci[i], label=f"Row={i}", color=COLORS[i], alpha=0.8)
        plt.title(f"{rootname}- all rows that are not in PSA/WCA", size=25)
        plt.ylim(-0.5, 4)
        plt.legend(loc="upper right")
        figname = os.path.join(outdir, f"{rootname}_{segment}_rows.png")
        plt.savefig(figname, bbox_inches="tight")
        print(f"Saved non-PSA/WCA rows: {figname}")
        plt.clf()
        val_c = c_stat(combined_dark, binned_sci, excluded_rows)

# TO DO, always hardcode this?
# Compare hardcoded vs variables
        #x0 = [0.007, 0.005]
        x0 = [0.1 * np.mean(binned_sci) / np.mean(lo_dark), 0.1 * np.mean(binned_sci) / np.mean(hi_dark)]
#        x0 = [lo_coeff, hi_coeff]
        res = minimize(fun_opt, x0, method="Nelder-Mead", tol=1e-6, 
                       args=([lo_dark, hi_dark], binned_sci, excluded_rows),
                       #bounds=[(1.e-8, None), (1.e-8, None)], options={'maxiter': 1000})
                       bounds=[(1.e-10, None), (1.e-10, None)], options={'maxiter': 1000})
        combined_dark1 = linear_combination([lo_dark, hi_dark], res.x)
        print("!!!!")
        print(res.x, lo_af["total_exptime"], hi_af["total_exptime"], sci_exp) 
        sh = np.shape(combined_dark1)
        plt.imshow(combined_dark1, aspect="auto", origin="lower",
                   extent=[0, sh[1], 0, sh[0]])
        plt.colorbar(label="Counts/s")
        plt.title(f"{rootname} Combined Superdark", size=25)
        figname = os.path.join(outdir, f"{rootname}_{segment}_combined_dark.png")
        plt.savefig(figname, bbox_inches="tight")
        print(f"Saved combined dark image: {figname}")
        plt.clf()

        # Use just one row outside of extraction regions
        ncols = 3
        nrows = int(np.ceil(len(bkg_rows)/ncols))
        fig, axes0 = plt.subplots(nrows, ncols, figsize=(25,15))
        axes = axes0.flatten()
        for i in range(len(bkg_rows)):
            row = bkg_rows[i]
            ax = axes[i]
            smoothy, smoothx = smooth_array(binned_sci[row], 25)
            avg = np.average(binned_sci[row])
            ax.plot(binned_sci[row], color="lightgrey", 
                    label=f"Sci row", alpha=0.8)
            ax.plot(combined_dark1[row], color=COLORS[0], 
                    label=f"Predicted Bkgd", alpha=0.8)
            ax.plot(smoothx, smoothy, color=COLORS[1], label="Smoothed Sci", alpha=0.8)
#            ax.axhline(avg, lw=2, color=COLORS[2], label=f"Avg Bkgd={avg:.1f}")
            ax.set_title(f"Row = {row}", size=15)
            ax.set_xlabel("X")
            ax.set_ylabel("Counts")
            ax.legend(loc="upper right")
        fig.suptitle("Predicted vs. Actual Dark across non PSA/WCA rows", size=25)
        figname = os.path.join(outdir, f"{rootname}_{segment}_predicted_dark_bkgd.png")
        plt.savefig(figname, bbox_inches="tight")
        print(f"Saved actual vs predicted darks: {figname}")
        plt.clf()
        
        # Use just one row outside of extraction regions
        ncols = 1
        nrows = int(np.ceil(len(apertures["PSA"])/ncols))
        fig, axes0 = plt.subplots(nrows, ncols, figsize=(25,15))
        axes = axes0.flatten()
        for i in range(len(apertures["PSA"])):
            row = apertures["PSA"][i]
            ax = axes[i]
            smoothy, smoothx = smooth_array(binned_sci[row], 25)
            avg = np.average(binned_sci[row])
            ax.plot(binned_sci[row], label="Sci row", color="lightgrey", alpha=0.8)
            ax.plot(combined_dark1[row], label=f"Predicted dark", color=COLORS[0], alpha=0.8)
            ax.plot(smoothx, smoothy, color=COLORS[1], label="Smoothed Sci", alpha=0.8)
#            ax.axhline(avg, lw=2, color=COLORS[2], label=f"Avg Sci={avg:.1f}")
            ax.set_title(f"Row = {row}", size=15)
            ax.set_xlabel("X")
            ax.set_ylabel("Counts")
            ax.legend(loc="upper right")
            ax.set_ylim(bottom=-1)
        figname = os.path.join(outdir, f"{rootname}_{segment}_predicted_dark_sci.png")
        fig.suptitle(f"{rootname}: Predicted dark and actual science across PSA rows", size=25)
        plt.savefig(figname, bbox_inches="tight")
        print(f"Saved predicted dark vs science: {figname}")
        plt.clf()

#       These two files below are used for sanity checks
        noise = copy.deepcopy(lo_af.tree)
        noise[pha_str] = combined_dark1[apertures["PSA"]]
        noise["scaling"] = res.x
        outfile = os.path.join(outdir, f"{rootname}_noise.asdf")
        noise_af = asdf.AsdfFile(noise)
        noise_af.write_to(outfile)
        print(f"Wrote {outfile}")

        signal = copy.deepcopy(lo_af.tree)
        signal[pha_str] = binned_sci[apertures["PSA"]]
        outfile = os.path.join(outdir, f"{rootname}_signal.asdf")
        signal_af = asdf.AsdfFile(signal)
        signal_af.write_to(outfile)
        print(f"Wrote {outfile}")

        noise_comp = copy.deepcopy(lo_af.tree)
        noise_comp[pha_str] = combined_dark1
        #combined_exptime = (lo_exptime * res.x[0]) + (hi_exptime * res.x[1])
        noise_comp["total_exptime"] = sci_exp
        outfile = os.path.join(outdir, f"{rootname}_noise_complete.asdf")
        noise_comp_af = asdf.AsdfFile(noise_comp)
        noise_comp_af.write_to(outfile)
        print(f"Wrote {outfile}")

    lo_af.close()
    hi_af.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datadir",
                        help="Path to science corrtags")
    parser.add_argument("--hv", default=None,
                        help="HV to filter corrtags by")
    parser.add_argument("--segment", default=None,
                        help="Segment to filter corrtags by")
    parser.add_argument("--lo", dest="lo_darkname",
                        help="Name of low activity superdark")
    parser.add_argument("--hi", dest="hi_darkname",
                        help="Name of high activity superdark")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Name of output directory")
    parser.add_argument("--binned", default=False,
                        action="store_true",
                        help="Toggle to indicate that supplied superdarks are binned")
    args = parser.parse_args()
    corrtags = glob.glob(os.path.join(args.datadir, "*corrtag*fits"))

    predict_dark(corrtags, args.lo_darkname, args.hi_darkname, args.segment, args.hv, args.outdir, args.binned)