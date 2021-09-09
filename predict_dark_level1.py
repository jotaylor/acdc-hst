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
# TO DO, hardcoded one pha range
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


def c_stat(combined_dark, binned_sci, excluded_rows=[2,3,4,7,8]):
    csum = 0.
    for y in range(combined_dark.shape[0]):
        for x in range(combined_dark.shape[1]):
            if binned_sci[y,x] > 0 and y not in excluded_rows and combined_dark[y,x] > 0:
                csum = csum + 2. * (combined_dark[y,x] - binned_sci[y,x] + binned_sci[y,x] * (np.log(binned_sci[y,x]) - np.log(combined_dark[y,x])))
            elif y not in excluded_rows:
                csum += np.abs(combined_dark[y,x]) * 2.

    return csum


def main(corrtags, lo_darkname, hi_darkname, segment=None, hv=None, outdir="."):
    Lo = Superdark.from_asdf(lo_darkname)
    Hi = Superdark.from_asdf(hi_darkname)
    lo_binnedname = lo_darkname.replace(".asdf", "_phabinned.asdf")
    hi_binnedname = hi_darkname.replace(".asdf", "_phabinned.asdf")
    Lo.bin_superdark(RESEL[0], RESEL[1], pha_bins=PHA_INCL_EXCL, outfile=lo_binnedname)
    Hi.bin_superdark(RESEL[0], RESEL[1], pha_bins=PHA_INCL_EXCL, outfile=hi_binnedname)
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
    for item in corrtags:
        if segment is not None:
            file_segment = fits.getval(item, "segment")
            if hv is not None:
                file_hv = fits.getval(item, f"HVLEVEL{file_segment[-1]}", 1)
                if file_segment != segment.upper() and str(file_hv) != str(hv):
                    print(f"File does not match required settings: {item}")
                    continue
        rootname0 = fits.getval(item, "rootname")
        rootname = rootname0.lower()
        binned_sci, nevents = bin_science(item, binning)
        plt.imshow(binned_sci, aspect="auto", origin="lower")
        plt.savefig("binned_sci.png")
        plt.clf()
        sci_exp = fits.getval(item, "exptime", 1)
        # Initial guess is 0.5 contribution for each superdark
        lo_coeff = 0.5 / (lo_af["total_exptime"] / sci_exp)
        hi_coeff = 0.5 / (hi_af["total_exptime"] / sci_exp)
        combined_dark = linear_combination([lo_dark, hi_dark], [lo_coeff, hi_coeff])
# this is science extraction regions (use xtractab) for sci exp. LP
# also exclude wavecals
        excluded_rows = [2,3,4,7,8]
        for i in range(binned_sci.shape[0]):
            if i not in excluded_rows:
                plt.plot(binned_sci[i], label=i)
        plt.legend()
        plt.savefig("rows.png")
        plt.clf()
        val_c = c_stat(combined_dark, binned_sci, excluded_rows)

# TO DO, always hardcode this?
# Compare hardcoded vs variables
        x0 = [0.007, 0.005]
#        x0 = [lo_coeff, hi_coeff]
        res = minimize(fun_opt, x0, method="Nelder-Mead", tol=1e-6, 
                       args=([lo_dark, hi_dark], binned_sci, excluded_rows))
        combined_dark1 = linear_combination([lo_dark, hi_dark], res.x)
        plt.imshow(combined_dark1, aspect="auto", origin="lower")
        plt.savefig("combined_dark.png")
        plt.clf()

        sci_row = 9
        plt.plot(binned_sci[sci_row], label="Outside of sci extr")
        plt.plot(combined_dark1[sci_row], label="Predicted dark level")
        plt.xlabel("X")
        plt.ylabel("Number of photons")
        plt.legend()
        plt.savefig("predicted_dark.png")
        plt.clf()

#       These two files below are used for sanity checks
        noise = copy.deepcopy(lo_af.tree)
        noise[pha_str] = combined_dark1[sci_row]
        outfile = os.path.join(outdir, f"{rootname}_noise.asdf")
        noise_af = asdf.AsdfFile(noise)
        noise_af.write_to(outfile)
        print(f"Wrote {outfile}")

        signal = copy.deepcopy(lo_af.tree)
        signal[pha_str] = binned_sci[sci_row]
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
    args = parser.parse_args()
    corrtags = glob.glob(os.path.join(args.datadir, "*corrtag*fits"))

    main(corrtags, args.lo_darkname, args.hi_darkname, args.segment, args.hv, args.outdir)
