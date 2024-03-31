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
import datetime
from collections import defaultdict
import itertools
import copy
import argparse
import os
import glob
import asdf
from scipy.optimize import minimize
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
dirname = os.path.dirname(__file__)
stylesheet = os.path.join(dirname, "analysis", "niceplot.mplstyle")
plt.style.use(stylesheet)
from matplotlib.backends.backend_pdf import PdfPages
import warnings
from astropy.units import UnitsWarning
warnings.filterwarnings("ignore", category=UnitsWarning)
import pandas as pd
pd.options.mode.chained_assignment = None

from acdc.superdark.cos_fuv_superdark import Superdark, get_gsag_holes
from acdc.database.calculate_dark import get_aperture_region
from acdc.analysis.compare_backgrounds import smooth_array
from acdc.utils.utils import bin_coords, timefunc, get_psa_wca


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
COLORS = ["#9970ab", "#5aae61", "#d95f02", "#e7298a", "#66a61e", "#a6cee3", 
          "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", 
          "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#a6cee3"]
# Hydrogen Lyman alpha wavelength limits
LYA_LO = 1213.67 # Angstroms
LYA_HI = 1217.67 # Angstroms
# Oxygen I wavelength limits
OI_LO = 1300.5
OI_HI = 1307.0

def measure_gsag(corrtag, row_threshold=.3):
    if fits.getval(corrtag, "life_adj") != 1:
        lp1_interpolate = False
    else:
        lp1_interpolate = None
    do_interpolate = True
    gsagtab = fits.getval(corrtag, "gsagtab")
    if "$" in gsagtab:
        lref = os.environ["lref"]
        gsagtab = os.path.join(lref, gsagtab.split("$")[1])
    segment = fits.getval(corrtag, "segment")
    hv = fits.getval(corrtag, f"HVLEVEL{segment[-1]}", 1)
    mjdend = fits.getval(corrtag, "EXPEND", 1)
    gsag_holes = get_gsag_holes(gsagtab, segment, hv)
    gsag_df = gsag_holes[gsag_holes["DATE"] < mjdend]
    if len(gsag_df) == 0:
        return do_interpolate, lp1_interpolate
    def prod(x0, y0, dx, dy):
        out = list(itertools.product(range(x0, x0+dx+1), range(y0, y0+dy+1)))
        return out
    coords = []
    for i in range(len(gsag_df)):
        s = gsag_df.iloc[i]
        c = prod(s["X0"], s["Y0"], s["DX"], s["DY"])
        coords += c

    #gsag_df["COORDS"] = gsag_df.apply(lambda row: prod(row['X0'], row['Y0'], row['DX'], row['DY']), axis=1)
    #coords0 = gsag_df['COORDS'].to_numpy()
    #coords = []
    #for item in coords0: #each element is a list
    #    coords += item
#    x = [x[0] for x in coords]
#    y = [x[1] for x in coords]
#    plt.plot(x, y, "ko")
#    plt.show()
    coords = list(set(coords))
    coords_x = np.array([x[0] for x in coords])
    coords_y = np.array([x[1] for x in coords])
    # Active Area limits from BRFTAB x1u1459il_brf.fits (essentially static)
    aa_lims = {"FUVA": (1060, 296, 15250, 734), #xstart, ystart, xend, yend
               "FUVB": (809,  360, 15182, 785)} #xstart, ystart, xend, yend
    seg_aa_lims = aa_lims[segment]
    pix_row = seg_aa_lims[2] - seg_aa_lims[0]
    for y in range(seg_aa_lims[1], seg_aa_lims[3]+1):
        inds = np.where((coords_y == y) & (coords_x > seg_aa_lims[0]) & (coords_x < seg_aa_lims[2]))
        frac = len(inds[0])/pix_row
        if frac >= row_threshold:
            do_interpolate = False
            break

    print(f"Interpolate={do_interpolate} for {corrtag}")
    return do_interpolate, lp1_interpolate

#    ngsag = sum(gsag_df["DX"] * gsag_df["DY"])
#    npix = 16777216 #16384*1024
#    det_threshold = (row_threshold*16384)/npix
#    gsag_frac = ngsag/npix
#    do_interpolate = True
#    if gsag_frac >= det_threshold:
#        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#        print(corrtag)
#        print(f"Number of gain-sagged pixels exceeds threshold!")
#        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#        do_interpolate = False

#    data = fits.getdata(corrtag)
#    gsaginds = np.where(data["dq"]&8192 == 8192)
#    x = data["xfull"][gsaginds]
#    y = data["yfull"][gsaginds]
#    x = x.astype(int)
#    y = y.astype(int)
#    coords = [(x[i], y[i]) for i in range(len(x))]
#    uniq_coords = list(set(coords))
#    ngsag = len(uniq_coords)
#    nevents = len(data["xfull"])
#    print((ngsag/nevents)*100.)
#    npix = 16777216 #16384*1024
#    det_threshold = (row_threshold*16384)/npix
#    if ngsag/npix >= det_threshold:
#        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
#        print(f"Number of gain-sagged pixels exceeds threshold!")
#        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")

#@timefunc
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
    # Uncomment the lines below for debugging on the convergence process
#    print(coeffs)
#    print(cval)
#    print("")
    
    return cval


#@timefunc
def linear_combination(darks, coeffs):
    """Scale and combine superdarks.

    Scale multiple superdarks (typically 1 quiescent and 1 active superdark),
    then add them together to create a scaled model superdark for a particular
    science exposure.

    Arguments:
        darks (list): List of superdark images to scale and combine.
        coeffs (list): List of scale factors to apply to each superdark.
    Returns:
        combined (array): Scaled and combined model superdark.
    """

    combined = darks[0] * coeffs[0]
    for i in range(1, len(darks)):
        combined += darks[i] * coeffs[i]
    return combined 


def check_superdarks(binned_superdarks):
    """Ensure all input superdarks were binned with identical bin sizes.

    Arguments: 
        binned_superdarks (list or array-like): List of binned superdarks.
    """

    bad = False
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
            "phastart", "phaend"]
    dark_vals = defaultdict(list)
    for darkfile in binned_superdarks:
        dark_af = asdf.open(darkfile)
        for k in keys:
            dark_vals[k].append(dark_af[k])
        dark_af.close()
    for k in keys:
        if len(set(dark_vals[k])) != 1:
            print(f"WARNING!!!! Key {k} does not match for all superdarks")
            bad = True
    if bad is True:
        for darkfile in binned_superdarks:
            print(darkfile)
        raise KeyError("Not all input superdarks are binned identically")


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

def bin_science(corrtag, b, fact=1, exclude_airglow=False):
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
    segment = fits.getval(corrtag, "segment")
    cenwave = fits.getval(corrtag, "cenwave")
    lp = fits.getval(corrtag, "life_adj")
    aperture_regions = get_aperture_region(cenwave=cenwave, aperture=aperture, segments=[segment], life_adj=[lp])
    box = aperture_regions[segment][f"lp{lp}_{aperture.lower()}_{cenwave}"]
    xmin0, xmax0, ymin0, ymax0 = box
    data = fits.getdata(corrtag)
    inds0 = np.where(
                     (data["pha"] >= b["phastart"]) & 
                     (data["pha"] < b["phaend"]) &
                     (data["xcorr"] >= b["xstart"]) &
                     (data["xcorr"] < b["xend"]) & 
                     (data["ycorr"] >= b["ystart"]) &
                     (data["ycorr"] < b["yend"]))
    filtered0 = data[inds0]
    if exclude_airglow is True:
        print(f"Excluding pixels affected by LyAlpha & OI airglow") 
        bad = (filtered0["wavelength"] > LYA_LO) &\
              (filtered0["wavelength"] < LYA_HI) &\
              (filtered0["wavelength"] > OI_LO) &\
              (filtered0["wavelength"] < OI_HI)
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
    
    xbinned, ybinned = bin_coords(filtered["xcorr"], filtered["ycorr"],
                            b["bin_x"], b["bin_y"], b["xstart"], b["ystart"],
                            make_int=True)
    
    for i in range(len(xbinned)):
        x = xbinned[i]
        y = ybinned[i]
        binned[y,x] += filtered["epsilon"][i]
        nevents[y,x] += 1

    return binned, nevents


#@timefunc
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
            elif y not in excluded_rows and combined_dark[y,x] > 0:
                csum += np.abs(combined_dark[y,x]) * 2.
            elif y not in excluded_rows and combined_dark[y,x] <= 0:
                csum += np.inf

    return csum


def check_existing(filename, overwrite=False):
    if os.path.exists(filename) and overwrite is False:
        print(f"WARNING: File {filename} already exists and overwrite is False, skipping...")
        return True
    else:
        return False


def predict_dark(corrtags, superdarks, segment=None, hv=None, 
                 outdir=".", binned=False, overwrite=False, exclude_airglow=False,
                 x_bin=RESEL[0]*2, y_bin=RESEL[1]*2):
    """Model the superdark for each input science corrtag.

    Use the active and quiescent superdarks of the
    appropriate segment+HV combination to determine the best predicted
    superdark model for each science corrtag.

    Args:
        corrtags (list or array-like): Predict the model superdark for these corrtags.
        superdarks (list or array-like): Superdarks to use for modeling. 
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
    
    # Ensure all input superdarks are binned identically 
    check_superdarks(superdarks)

    for item in corrtags:
        do_interpolate_gsag, do_interpolate_lp1 = measure_gsag(item)
        if do_interpolate_gsag is True:
            break

    # If superdarks are not yet binned, bin them.
    if binned is False:
        binned_superdarks = []
        for darkfile in superdarks:
            S = Superdark.from_asdf(darkfile, overwrite=overwrite)
            binnedname0 = darkfile.replace(".asdf", "_binned.asdf")
            binnedname = os.path.join(outdir, binnedname0)
            binned_superdarks.append(binnedname)
            S.screen_counts(verbose=False, interpolate=True)
            S.bin_superdark(RESEL[0]*2, RESEL[1]*2, pha_bins=PHA_INCL_EXCL,
                             writefile=False)
            S.screen_counts(verbose=False, sigma=5, interpolate=True)
            if S.fixed_gsag is False and do_interpolate_gsag is True:
                S.fix_gsag(lp1_interpolate=do_interpolate_lp1)
            S.write_superdark(user_outfile=binnedname)
            S.plot_superdarks()
    else:
        binned_superdarks = superdarks
        for darkfile in superdarks:
            S = Superdark.from_asdf(darkfile, overwrite=True)
            S.screen_counts(verbose=False, sigma=5, interpolate=True)
            if S.fixed_gsag is False and do_interpolate_gsag is True:
                S.fix_gsag(lp1_interpolate=do_interpolate_lp1)
            S.write_superdark(user_outfile=darkfile)
            S.plot_superdarks()

    # Plot the binned superdarks 
    pdfname = os.path.join(outdir, "all_binned_superdarks.pdf")
    pdf = PdfPages(pdfname)
    dark_ims = []
    superdark_exptimes = []
    for darkfile in binned_superdarks:
        dark_af = asdf.open(darkfile)
        superdark_exptimes.append(dark_af["total_exptime"])
        binning = get_binning_pars(dark_af)
        pha_str = f"pha{PHA_INCL_EXCL[0]}-{PHA_INCL_EXCL[1]}"
        dark = dark_af[pha_str]
        this = type(dark.data)
        dark_im = copy.deepcopy(dark)
        dark_ims.append(dark_im)
        dark_hv = dark_af["hv"]
        dark_segment = dark_af["segment"]
        fig,ax = plt.subplots(1, 1, figsize=(20,8))
        im = ax.imshow(dark, aspect="auto", origin="lower", interpolation="nearest")
        fig.colorbar(im, label="Counts", pad=0.01)
        ax.set_title(f"{os.path.basename(darkfile)}\n{dark_segment} {dark_hv} superdark", size=25)
        pdf.savefig(fig, bbox_inches="tight")
        dark_af.close()
    pdf.close()
    print(f"Wrote {pdfname}")

    # Now go through each input science corrtag
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
        excluded_rows, apertures = get_psa_wca(segment, cenwave, lp, binning)
        rootname0 = fits.getval(item, "rootname")
        rootname = rootname0.lower()
        binned_sci, nevents = bin_science(item, binning, segment, cenwave, lp, exclude_airglow=exclude_airglow)
#        vmin = np.mean(binned_sci) - 1*np.std(binned_sci)
#        if vmin < 0:
#            vmin = 0
        vmin = 0
        vmax = np.median(binned_sci) + np.median(binned_sci)*2.
        if vmax < 6:
            vmax = 6
        sh = np.shape(binned_sci)
        im = ax.imshow(binned_sci, aspect="auto", origin="lower", 
                vmin=vmin, vmax=vmax, extent=[0, sh[1], 0, sh[0]], interpolation="nearest")
        ax.axhline(apertures["PSA"][0], lw=2, color=COLORS[3], label="PSA")
        ax.axhline(apertures["PSA"][-1]+1, lw=2, color=COLORS[3])
        ax.axhline(apertures["WCA"][0], lw=2, ls="dashed", color=COLORS[3], label="WCA")
        ax.axhline(apertures["WCA"][-1]+1, lw=2, ls="dashed", color=COLORS[3])
        ax.legend(loc="upper right")
        fig.colorbar(im, label="Counts", pad=0.01)
        figname = os.path.join(outdir, f"{rootname}_{segment}_binned_sci.png")
        ax.set_title(f"{rootname} {segment} Binned Science Image", size=25)
        ax.set_xlabel("X (binned)")
        ax.set_ylabel("Y (binned)")
        exists = check_existing(figname, overwrite)
        if not exists:
            fig.savefig(figname, bbox_inches="tight")
            print(f"Saved binned science image: {figname}")
        plt.clf()
        sci_exp = fits.getval(item, "exptime", 1)

        nrows = np.shape(binned_sci)[0]
        all_rows = np.arange(nrows)
        bkg_rows = list(set(all_rows) - set(excluded_rows))

        # This plot is not really that useful nowadays.
#        # Plot the counts in each row that are used for scaling-
#        # that is, all non-PSA and non-WCA rows.        
#        fig, ax = plt.subplots(1, 1, figsize=(20, 8))
#        for i in range(binned_sci.shape[0]):
#            if i not in excluded_rows:
#                ax.plot(binned_sci[i], label=f"Row={i}", color=COLORS[i], alpha=0.8)
#        ax.set_title(f"{rootname}- all rows that are not in PSA/WCA", size=25)
#        ax.set_ylim(-0.5, 4)
#        ax.legend(loc="upper right")
#        figname = os.path.join(outdir, f"{rootname}_{segment}_rows.png")
#        exists = check_existing(figname, overwrite)
#        if not exists:
#            fig.savefig(figname, bbox_inches="tight")
#            print(f"Saved non-PSA/WCA rows: {figname}")
#        plt.clf()

        # Make an initial guess at the combined superdark which is just an equal
        # contribution from each binned superdark.
        equal = 1 / len(binned_superdarks)
        superdark_coeffs = [equal / (exptime / sci_exp) for exptime in superdark_exptimes]
        combined_dark = linear_combination(dark_ims, superdark_coeffs)

        # Now actually determine the best combination of all superdarks 
# TO DO, always hardcode this?
# Compare hardcoded vs variables
        #x0 = [0.007, 0.005]
        x0 = [0.1 * np.mean(binned_sci) / np.mean(dark_im) for dark_im in dark_ims]
#        x0 = superdark_coeffs
        res = minimize(fun_opt, x0, method="Nelder-Mead", tol=1e-6, 
                       args=(dark_ims, binned_sci, excluded_rows),
                       #bounds=[(1.e-8, None), (1.e-8, None)], options={'maxiter': 1000})
                       bounds=[(1.e-10, None) for x in binned_superdarks], options={'maxiter': 500})
        combined_dark1 = linear_combination(dark_ims, res.x)
        #print(res.x, lo_af["total_exptime"], hi_af["total_exptime"], sci_exp)
        now2 = datetime.datetime.now()
        print(now2)

        # Plot the model superdark
        sh = np.shape(combined_dark1)
        vmin = np.median(combined_dark1) - np.median(combined_dark1)*0.5
        if vmin < 0:
            vmin = 0
        vmax = np.median(combined_dark1) + np.median(combined_dark1)*2.
        fig, ax = plt.subplots(1, 1, figsize=(20, 8))
        im = ax.imshow(combined_dark1, aspect="auto", origin="lower",
                   extent=[0, sh[1], 0, sh[0]], interpolation="nearest",
                   vmin=vmin, vmax=vmax)
        fig.colorbar(im, label="Counts", pad=0.01)
        ax.set_title(f"{rootname} Best Model Superdark", size=25)
        ax.set_xlabel("X (binned)")
        ax.set_ylabel("Y (binned)")
        figname = os.path.join(outdir, f"{rootname}_{segment}_combined_dark.png")
        exists = check_existing(figname, overwrite)
        if not exists:
            fig.savefig(figname, bbox_inches="tight")
            print(f"Saved combined dark image: {figname}")
        plt.clf()

        # Plot the predicted vs. actual dark rates for all non-PSA and non-WCA rows
        ncols = 3
        nrows = int(np.ceil(len(bkg_rows)/ncols))
        fig, axes0 = plt.subplots(nrows, ncols, figsize=(25, nrows*6))
        fig.subplots_adjust(hspace=0.3, wspace=0.15)
        axes = axes0.flatten()
        for i in range(len(bkg_rows)):
            row = bkg_rows[i]
            ax = axes[i]
            smoothy, smoothx = smooth_array(binned_sci[row], 25)
            avg = np.average(binned_sci[row])
            avg_smoothy = np.average(smoothy)
            diff = np.abs(avg_smoothy - np.average(combined_dark1[row]))
            if diff > 2*avg_smoothy:
                ax.annotate("WARNING!\nBAD FIT!", (.05, .85), color=COLORS[2], xycoords="axes fraction")    
            ax.plot(binned_sci[row], color="lightgrey", 
                    label=f"Science", alpha=0.8)
            ax.plot(combined_dark1[row], color=COLORS[0], 
                    label=f"Predicted Bkgd.", alpha=0.8)
            ax.plot(smoothx, smoothy, color=COLORS[1], label="Smoothed Sci.", alpha=0.8)
            smoothy_max = np.max(smoothy) * 2.
            ylim = ax.get_ylim()
            if ylim[1] > smoothy_max:
                ax.set_ylim(top=smoothy_max)
#            ax.axhline(avg, lw=2, color=COLORS[2], label=f"Avg Bkgd={avg:.1f}")
            ax.set_title(f"Binned row = {row}", size=18)
            ax.set_xlabel("X (binned)")
            ax.set_ylabel("Counts")
        axes[2].legend(loc="upper right", fontsize=18) 
        #handles, labels = ax.get_legend_handles_labels()
        #fig.legend(handles, labels, bbox_to_anchor=(.9, .91), loc="lower right", fontsize=18) 
        fig.suptitle(f"Predicted vs. actual dark for non-PSA/WCA rows\n{rootname} {segment}", size=25, y=.95)
        figname = os.path.join(outdir, f"{rootname}_{segment}_predicted_dark_bkgd.png")
        exists = check_existing(figname, overwrite)
        if not exists:
            plt.savefig(figname, bbox_inches="tight")
            print(f"Saved actual vs predicted darks: {figname}")
        plt.clf()
        
        # Plot the predicted vs. actual dark rates for all PSA rows
        ncols = 1
        nrows = int(np.ceil(len(apertures["PSA"])/ncols))
        fig, axes0 = plt.subplots(nrows, ncols, figsize=(22, nrows*6))
        fig.subplots_adjust(hspace=0.3)
        axes = axes0.flatten()
        for i in range(len(apertures["PSA"])):
            row = apertures["PSA"][i]
            ax = axes[i]
            smoothy, smoothx = smooth_array(binned_sci[row], 25)
            avg = np.average(binned_sci[row])
            ax.plot(binned_sci[row], label="Science", color="lightgrey", alpha=0.8)
            ax.plot(combined_dark1[row], label=f"Predicted Bkgd.", color=COLORS[0], alpha=0.8)
            ax.plot(smoothx, smoothy, color=COLORS[1], label="Smoothed Sci.", alpha=0.8)
#            ax.axhline(avg, lw=2, color=COLORS[2], label=f"Avg Sci={avg:.1f}")
            ax.set_title(f"Binned row = {row}", size=20)
            ax.set_xlabel("X (binned)")
            ax.set_ylabel("Counts")
            ax.set_ylim(bottom=-.75)
            smoothy_max = np.max(smoothy) * 2.
            ylim = ax.get_ylim()
            if ylim[1] > smoothy_max:
                ax.set_ylim(top=smoothy_max)
        axes[0].legend(loc="upper right", fontsize=18) 
        #handles, labels = ax.get_legend_handles_labels()
        #fig.legend(handles, labels, bbox_to_anchor=(.9, .91), loc="lower right", fontsize=18) 
        figname = os.path.join(outdir, f"{rootname}_{segment}_predicted_dark_sci.png")
        fig.suptitle(f"Predicted vs. actual dark across PSA rows\n{rootname} {segment}", size=25, y=.97)
        exists = check_existing(figname, overwrite)
        if not exists:
            plt.savefig(figname, bbox_inches="tight")
            print(f"Saved predicted dark vs science: {figname}")
        plt.clf()

#       These two files below are used for sanity checks
#        noise = copy.deepcopy(lo_af.tree)
#        noise[pha_str] = combined_dark1[apertures["PSA"]]
#        noise["scaling"] = res.x
#        outfile = os.path.join(outdir, f"{rootname}_{segment}_noise.asdf")
#        noise_af = asdf.AsdfFile(noise)
#        exists = check_existing(outfile, overwrite)
#        if not exists:
#            noise_af.write_to(outfile)
#            print(f"Wrote {outfile}")
#
#        signal = copy.deepcopy(lo_af.tree)
#        signal[pha_str] = binned_sci[apertures["PSA"]]
#        outfile = os.path.join(outdir, f"{rootname}_{segment}_signal.asdf")
#        signal_af = asdf.AsdfFile(signal)
#        exists = check_existing(outfile, overwrite)
#        if not exists:
#            signal_af.write_to(outfile)
#            print(f"Wrote {outfile}")
        
        dark_af = asdf.open(binned_superdarks[0])
        noise_comp = copy.deepcopy(dark_af.tree)
        noise_comp[pha_str] = combined_dark1
        #combined_exptime = (lo_exptime * res.x[0]) + (hi_exptime * res.x[1])
        noise_comp["total_exptime"] = sci_exp
        outfile = os.path.join(outdir, f"{rootname}_{segment}_noise_complete.asdf")
        noise_comp_af = asdf.AsdfFile(noise_comp)
        exists = check_existing(outfile, overwrite)
        if not exists:
            noise_comp_af.write_to(outfile)
            print(f"Wrote {outfile}")
        dark_af.close()

