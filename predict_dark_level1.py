
import copy
import argparse
import os
import glob
import asdf
from scipy.optimize import minimize
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

from make_clean_superdark import bin_corrtag
from cos_fuv_superdark import Superdark
from calculate_dark import get_1291_box

RESEL = [6, 10]
PHA_INCLUSIVE = [2, 23]
PHA_INCL_EXCL = [2, 24]

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
    bad = False
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
            "phastart", "phaend"]
    for k in keys:
        if af1[k] != af2[k]:
            print(f"WARNING!!!! Key {k} does not match for both superdarks")
            bad = True
    assert bad == False, "Cannot continue until discrepancies are resolved"


def get_binning_pars(af):
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
            "phastart", "phaend"]
    binning = {}
    for k in keys:
        binning[k] = af[k]

    return binning

def bin_science(corrtag, b, fact=1):
    data = fits.getdata(corrtag)
    phainds = np.where((data["pha"] >= b["phastart"]) & 
                       (data["pha"] <= b["phaend"]))
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


def bin_coords(xs, ys, bin_x, bin_y, xstart=0, ystart=0):
    xsnew = xs // bin_x
    ysnew = ys // bin_y
    bin_xstart = xstart // bin_x
    bin_ystart = ystart // bin_y
    xsnew -= bin_xstart
    ysnew -= bin_ystart
    return xsnew, ysnew

def get_excluded_rows(segment, lp, binning):
    psa = 0
    excluded_rows = np.array(())
    for aperture in ["PSA", "WCA"]:
        psa_1291 = get_1291_box(aperture=aperture)
        box = psa_1291[segment][f"lp{lp}_psa_1291"]
        xmin0, xmax0, ymin0, ymax0 = box
        xsnew, ysnew = bin_coords(np.array([xmin0, xmax0]), np.array([ymin0, ymax0]), binning["bin_x"], binning["bin_y"], binning["xstart"], binning["ystart"])
        psa_xmin, psa_xmax = xsnew
        psa_ymin, psa_ymax = ysnew
        rows = np.arange(psa_ymin, psa_ymax+1)
        if aperture == "PSA":
            psa = rows
        excluded_rows = np.concatenate((excluded_rows, rows))
    excluded_rows = excluded_rows.astype(int)
    psa = psa.astype(int) 
    
    return excluded_rows, psa

def c_stat(combined_dark, binned_sci, excluded_rows):
    csum = 0.
    for y in range(combined_dark.shape[0]):
        for x in range(combined_dark.shape[1]):
            if binned_sci[y,x] > 0 and y not in excluded_rows and combined_dark[y,x] > 0:
                csum = csum + 2. * (combined_dark[y,x] - binned_sci[y,x] + binned_sci[y,x] * (np.log(binned_sci[y,x]) - np.log(combined_dark[y,x])))
            elif y not in excluded_rows:
                csum += np.abs(combined_dark[y,x]) * 2.

    return csum


def main(corrtags, lo_darkname, hi_darkname, segment=None, hv=None, outdir=".", binned=False):
    if binned is False:
        Lo = Superdark.from_asdf(lo_darkname)
        Hi = Superdark.from_asdf(hi_darkname)
        lo_binnedname = lo_darkname.replace(".asdf", "_phabinned.asdf")
        hi_binnedname = hi_darkname.replace(".asdf", "_phabinned.asdf")
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
    plt.imshow(lo_dark, aspect="auto", origin="lower")
    plt.savefig("lo_dark.png")
    plt.clf()
    plt.imshow(hi_dark, aspect="auto", origin="lower")
    plt.savefig("hi_dark.png")
    plt.clf()
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

        lp = fits.getval(item, "life_adj")
        excluded_rows, psa = get_excluded_rows(segment, lp, binning)
        rootname0 = fits.getval(item, "rootname")
        rootname = rootname0.lower()
        binned_sci, nevents = bin_science(item, binning)
        plt.imshow(binned_sci, aspect="auto", origin="lower")
        plt.savefig(f"{rootname}_binned_sci.png")
        plt.clf()
        sci_exp = fits.getval(item, "exptime", 1)
        # Initial guess is 0.5 contribution for each superdark
        lo_coeff = 0.5 / (lo_af["total_exptime"] / sci_exp)
        hi_coeff = 0.5 / (hi_af["total_exptime"] / sci_exp)
        combined_dark = linear_combination([lo_dark, hi_dark], [lo_coeff, hi_coeff])
# this is science extraction regions (use xtractab) for sci exp. LP
# also exclude wavecals
        for i in range(binned_sci.shape[0]):
            if i not in excluded_rows:
                plt.plot(binned_sci[i], label=i)
        plt.legend()
        plt.savefig(f"{rootname}_rows.png")
        plt.clf()
        val_c = c_stat(combined_dark, binned_sci, excluded_rows)

# TO DO, always hardcode this?
# Compare hardcoded vs variables
        x0 = [0.007, 0.005]
#        x0 = [lo_coeff, hi_coeff]
        res = minimize(fun_opt, x0, method="Nelder-Mead", tol=1e-6, 
                       args=([lo_dark, hi_dark], binned_sci, excluded_rows),
                       bounds=[(0.0, None), (0, None)], options={'maxiter': 1000})
        combined_dark1 = linear_combination([lo_dark, hi_dark], res.x)
        plt.imshow(combined_dark1, aspect="auto", origin="lower")
        plt.savefig(f"{rootname}_combined_dark.png")
        plt.clf()

        # Use just one row outside of extraction regions
        row = psa[0]
        plt.plot(binned_sci[row], label=f"Binned science row={row}")
        plt.plot(combined_dark1[row], label=f"Predicted dark row={row}")
        plt.xlabel("X")
        plt.ylabel("Number of photons")
        plt.legend()
        plt.savefig(f"{rootname}_predicted_dark.png")
        plt.clf()

#       These two files below are used for sanity checks
        noise = copy.deepcopy(lo_af.tree)
        noise[pha_str] = combined_dark1[psa]
        outfile = os.path.join(outdir, f"{rootname}_noise.asdf")
        noise_af = asdf.AsdfFile(noise)
        noise_af.write_to(outfile)
        print(f"Wrote {outfile}")

        signal = copy.deepcopy(lo_af.tree)
        signal[pha_str] = binned_sci[psa]
        outfile = os.path.join(outdir, f"{rootname}_signal.asdf")
        signal_af = asdf.AsdfFile(signal)
        signal_af.write_to(outfile)
        print(f"Wrote {outfile}")

        noise_comp = copy.deepcopy(lo_af.tree)
        noise_comp[pha_str] = combined_dark1
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

    main(corrtags, args.lo_darkname, args.hi_darkname, args.segment, args.hv, args.outdir, args.binned)
