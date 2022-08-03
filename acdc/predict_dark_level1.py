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


# Hardcoded constant values
RESEL = [6, 10] # FUV resolution element size
# These are the PHAs that CalCOS uses for calibration, where the 
# lower and upper bounds are both inclusive 
PHA_INCLUSIVE = [2, 23] 
# These are the PHAs that CalCOS uses for calibration, where the 
# lower bound is inclusive and upper bounds is exclusive. This is
# the logic python uses!!
PHA_INCL_EXCL = [2, 24]
# Colorblind-safe palette below
COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6cee3", 
          "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", 
          "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#a6cee3"]
# Lyman alpha wavelength limits
LYA_LO = 1213.67 # Angstroms
LYA_HI = 1217.67 # Angstroms

def fun_opt(coeffs, darks, binned_sci, excluded_rows):
    """TODO- fill in

    Arguments:
        coeffs (list): The output coefficients.
        darks (list): A list of all superdarks to model, typicall [quiescent_dark, active_dark].
        binned_sci (array): Binned science image.
        excluded_rows (array): List of rows to exclude from modeling. Correspond
            to the PSA and WCA regions.
    Returns:
        cval (list): Coefficients? TODO
    """

    combined_dark = linear_combination(darks, coeffs)
    cval = c_stat(combined_dark, binned_sci, excluded_rows)
    
    return cval


def linear_combination(darks, coeffs):
    """Scale and combine superdarks.

    Scale multiple superdarks (typically 1 quiescent and 1 active superdark),
    then add them together to create a scaled model superdark for a particular
    science exposure.

    Arguments:
        darks (list): List of superdarks to scale and combine.
        coeffs (list): List of scale factors to apply to each superdark.
    Returns:
        combined (array): Scaled and combined model superdark.
    """

    combined = darks[0] * coeffs[0]
    for i in range(1, len(darks)):
        combined += darks[i] * coeffs[i]
    return combined 


def check_superdarks(af1, af2):
    """Ensure all input superdarks were binned with identical bin sizes.

    Arguments: 
        af1 (AsdfFile): ASDF file 1 to compare.
        af2 (AsdfFile): ASDF file 2 to compare.
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
    """Return spatial and PHA binning information for a superdark. 
    
    Arguments:
        af (AsdfFile): ASDF file to extract binning information from.
    Returns:
        binning (dict): Dictionary that describes the binning information
            in both the spatial and PHA dimensions.
    """
    
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
            "phastart", "phaend"]
    binning = {}
    for k in keys:
        binning[k] = af[k]

    return binning

def bin_science(corrtag, b, segment, cenwave, lp, fact=1, exclude_lya=False):
    """
    Given a corrtag with lists of events as a function of X, Y, and PHA,
    bin the data into an image using the same bin sizes as the superdark.

    Arguments:
        corrtag (str): Corrtag to bin.
        b (dict): Dictionary that describes the binning information
            in both the spatial and PHA dimensions. (Only spatial is used)
    Returns:
        binned (array): Binned science image.
        nevents (array): Number of events in each binned pixel.
    """
    
    aperture = "PSA"
    aperture_regions = get_aperture_region(cenwave=cenwave, aperture=aperture, segments=[segment], life_adj=[lp])
    box = aperture_regions[segment][f"lp{lp}_{aperture.lower()}_{cenwave}"]
    xmin0, xmax0, ymin0, ymax0 = box
    data = fits.getdata(corrtag)
    inds0 = np.where(
                     (data["pha"] >= b["phastart"]) & 
                     (data["pha"] < b["phaend"]) &
                     (data["xcorr"] > b["xstart"]) &
                     (data["xcorr"] < b["xend"]) & 
                     (data["ycorr"] > b["ystart"]) &
                     (data["ycorr"] < b["yend"]))
    filtered0 = data[inds0]
    if exclude_lya is True:
        print(f"Excluding pixels affected by Lyman Alpha airglow") 
        bad = (filtered0["wavelength"] > LYA_LO) &\
              (filtered0["wavelength"] < LYA_HI)
        #      (filtered0["ycorr"] > ymin0) &\
        #      (filtered0["ycorr"] < ymax0)
        inds1 = ~bad
        filtered = filtered0[inds1]
    else:
        filtered = filtered0

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
    igiven segment, cenwave, and lifetime position. Return the row indices in 
    the binned superdark coordinate system.

    Arguments:
        segment (str): Get PSA/WCA regions for this segment.
        cenwave (str): Get PSA/WCA regions for this cenwave.
        lp (str): Get PSA/WCA regions for this lifetime position.
        binning (dict): Dictionary that describes the binning information
            in both the spatial and PHA dimensions.
    Returns:
        excluded_rows (array): List of rows to exclude. Corresponds
            to the PSA and WCA regions.
        apertures (dict): The excluded rows for each aperture, PSA and WCA.
    """
    
    excluded_rows = np.array(())
    apertures = {"PSA": None, "WCA": None}
    for aperture in apertures:
        aperture_regions = get_aperture_region(cenwave=cenwave, aperture=aperture, segments=[segment])
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
    """TODO document

    Arguments:
        combined_ark (array): Scaled and combined model superdark.
        binned_sci (array): Binned science image.
        excluded_rows (array): List of rows to exclude. Corresponds
            to the PSA and WCA regions.
    
    Returns:
        csum (float): TODO
    """

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


def check_existing(filename, overwrite=False):
    if os.path.exists(filename) and overwrite is False:
        print(f"WARNING: File {filename} already exists and overwrite is False, skipping...")
        return True
    else:
        return False


def predict_dark(corrtags, lo_darkname, hi_darkname, segment=None, hv=None, 
                 outdir=".", binned=False, overwrite=False, exclude_lya=False):
    """Model the superdark for each input science corrtag.

    Use the active and quiescent superdarks of the
    appropriate segment+HV combination to determine the best predicted
    superdark model for each science corrtag.

    Args:
        corrtags (list or array-like): Predict the model superdark for these corrtags.
        lo_darkname (str): Quiescent superdark of the appropriate segment+HV combination.
        hi_darkname (str): Active superdark of the appropriate segment+HV combination.
        segment (str): (Optional) Process corrtags of this segment only.
        hv (str): (Optional) Process corrtags of this HV only.
        outdir (str): (Optional) Output directory to write products and plots.
            Default is current working directory.
        binned (Bool): (Optional) If True, indicates input superdarks are already
            binned. If False, input superdarks will be binned on the fly.
        overwrite (Bool): If True, overwrite any pre-existing products.
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if binned is False:
        Lo = Superdark.from_asdf(lo_darkname, overwrite=overwrite)
        Hi = Superdark.from_asdf(hi_darkname, overwrite=overwrite)
        lo_binnedname0 = lo_darkname.replace(".asdf", "_phabinned.asdf")
        hi_binnedname0 = hi_darkname.replace(".asdf", "_phabinned.asdf")
        lo_binnedname = os.path.join(outdir, lo_binnedname0)
        hi_binnedname = os.path.join(outdir, hi_binnedname0)
        Lo.screen_counts(verbose=False)
        Hi.screen_counts(verbose=False)
        Lo.bin_superdark(RESEL[0]*2, RESEL[1]*2, pha_bins=PHA_INCL_EXCL,
                         writefile=False)
        Hi.bin_superdark(RESEL[0]*2, RESEL[1]*2, pha_bins=PHA_INCL_EXCL,
                         writefile=False) 
        Lo.screen_counts(verbose=False, sigma=5, mask=True)
        Hi.screen_counts(verbose=False, sigma=5, mask=True)
        Lo.write_superdark
        Hi.write_superdark
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
    exists = check_existing(figname, overwrite)
    if not exists:
        fig.savefig(figname, bbox_inches="tight")
        print(f"Saved quiscent superdark figure: {figname}")
    fig,ax = plt.subplots(1, 1, figsize=(20,8))
    im = ax.imshow(hi_dark, aspect="auto", origin="lower")
    fig.colorbar(im, label="Counts")
    ax.set_title(f"{dark_segment} {dark_hv} active superdark", size=25)
    figname = os.path.join(outdir, f"{dark_segment}_{dark_hv}_hi_superdark.png")
    exists = check_existing(figname, overwrite)
    if not exists:
        fig.savefig(figname, bbox_inches="tight")
        print(f"Saved active superdark figure:{figname}")
    for item in corrtags:
        file_segment = fits.getval(item, "segment")
        file_hv = fits.getval(item, f"HVLEVEL{file_segment[-1]}", 1)
        if segment is not None:
            if file_segment != segment.upper() and str(file_hv) != str(hv):
                print(f"File does not match required settings: {item}")
                continue
        else:
            segment = file_segment
        if hv is not None:
            if file_hv != hv:
                print(f"File does not match required settings: {item}")
                continue

        cenwave = fits.getval(item, "cenwave")
        fig, ax = plt.subplots(1, 1, figsize=(20, 8))
        lp = fits.getval(item, "life_adj")
        excluded_rows, apertures = get_excluded_rows(segment, cenwave, lp, binning)
        rootname0 = fits.getval(item, "rootname")
        rootname = rootname0.lower()
        binned_sci, nevents = bin_science(item, binning, segment, cenwave, lp, exclude_lya=exclude_lya)
        print(binning)
#        vmin = np.mean(binned_sci) - 1*np.std(binned_sci)
#        if vmin < 0:
#            vmin = 0
#        vmax = np.mean(binned_sci) + 1*np.std(binned_sci)
        vmin = 0
        vmax = 10
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
        ax.set_title(f"{rootname} {segment} Binned Science Image", size=25)
        exists = check_existing(figname, overwrite)
        if not exists:
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
        exists = check_existing(figname, overwrite)
        if not exists:
            plt.savefig(figname, bbox_inches="tight")
            print(f"Saved non-PSA/WCA rows: {figname}")
        plt.clf()
#        val_c = c_stat(combined_dark, binned_sci, excluded_rows)

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
        #print(res.x, lo_af["total_exptime"], hi_af["total_exptime"], sci_exp) 
        sh = np.shape(combined_dark1)
        plt.imshow(combined_dark1, aspect="auto", origin="lower",
                   extent=[0, sh[1], 0, sh[0]])
        plt.colorbar(label="Counts/s")
        plt.title(f"{rootname} Combined Superdark", size=25)
        figname = os.path.join(outdir, f"{rootname}_{segment}_combined_dark.png")
        exists = check_existing(figname, overwrite)
        if not exists:
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
        fig.suptitle("Predicted vs. Actual Dark across non PSA/WCA rows\n{rootname} {segment}", size=25)
        figname = os.path.join(outdir, f"{rootname}_{segment}_predicted_dark_bkgd.png")
        exists = check_existing(figname, overwrite)
        if not exists:
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
        fig.suptitle("Predicted vs. Actual Dark across PSA rows\n{rootname} {segment}", size=25)
        exists = check_existing(figname, overwrite)
        if not exists:
            plt.savefig(figname, bbox_inches="tight")
            print(f"Saved predicted dark vs science: {figname}")
        plt.clf()

#       These two files below are used for sanity checks
        noise = copy.deepcopy(lo_af.tree)
        noise[pha_str] = combined_dark1[apertures["PSA"]]
        noise["scaling"] = res.x
        outfile = os.path.join(outdir, f"{rootname}_{segment}_noise.asdf")
        noise_af = asdf.AsdfFile(noise)
        exists = check_existing(outfile, overwrite)
        if not exists:
            noise_af.write_to(outfile)
            print(f"Wrote {outfile}")

        signal = copy.deepcopy(lo_af.tree)
        signal[pha_str] = binned_sci[apertures["PSA"]]
        outfile = os.path.join(outdir, f"{rootname}_{segment}_signal.asdf")
        signal_af = asdf.AsdfFile(signal)
        exists = check_existing(outfile, overwrite)
        if not exists:
            signal_af.write_to(outfile)
            print(f"Wrote {outfile}")

        noise_comp = copy.deepcopy(lo_af.tree)
        noise_comp[pha_str] = combined_dark1
        #combined_exptime = (lo_exptime * res.x[0]) + (hi_exptime * res.x[1])
        noise_comp["total_exptime"] = sci_exp
        outfile = os.path.join(outdir, f"{rootname}_{segment}_noise_complete.asdf")
        noise_comp_af = asdf.AsdfFile(noise_comp)
        exists = check_existing(outfile, overwrite)
        if not exists:
            noise_comp_af.write_to(outfile)
            print(f"Wrote {outfile}")

    lo_af.close()
    hi_af.close()

