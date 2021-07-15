import argparse
import os
import glob
import asdf
from scipy.optimize import minimize
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

from make_clean_superdark import bin_corrtag


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
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend",
            "phastart", "phaend"]
    bad = False
    binning = {}
    for k in keys:
        if af1[k] != af2[k]:
            print(f"WARNING!!!! Key {k} does not match for both superdarks")
            bad = True
        binning[k] = af1[k]
    assert bad == False, "Cannot continue until discrepancies are resolved"

    return binning


def bin_science(corrtag, b):
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
    binned = np.zeros((ydim, xdim))
    
    for i in range(len(filtered["xcorr"])):
        x = int((filtered["xcorr"][i] - b["xstart"]) // b["bin_x"])
        y = int((filtered["ycorr"][i] - b["ystart"]) // b["bin_y"])
        binned[y,x] += filtered["epsilon"][i]
    return binned


def c_stat(combined_dark, binned_sci, excluded_rows=[2,3,4,7,8]):
    csum = 0.
    for y in range(combined_dark.shape[0]):
        for x in range(combined_dark.shape[1]):
            if binned_sci[y,x] > 0 and y not in excluded_rows and combined_dark[y,x] > 0:
                csum = csum + 2. * (combined_dark[y,x] - binned_sci[y,x] + binned_sci[y,x] * (np.log(binned_sci[y,x]) - np.log(combined_dark[y,x])))
            elif y not in excluded_rows:
                csum += np.abs(combined_dark[y,x]) * 2.

    return csum


def main(corrtags, lo_darkname, hi_darkname):
    lo_af = asdf.open(lo_darkname)
    hi_af = asdf.open(hi_darkname)
    binning = check_superdarks(lo_af, hi_af)
# TO DO, fix hardcoded key
    lo_dark = lo_af["3-29"]
    hi_dark = hi_af["3-29"]
    plt.imshow(lo_dark, aspect="auto", origin="lower")
    plt.savefig("lo_dark.png")
    plt.clf()
    for item in corrtags:
        binned_sci = bin_science(item, binning)
        plt.imshow(binned_sci, aspect="auto", origin="lower")
        plt.savefig("binned_sci.png")
        plt.clf()
        sci_exp = fits.getval(item, "exptime", 1)
# TO DO, check with Andrei about this number
        lo_coeff = 0.5 / (lo_af["total_exptime"] / sci_exp)
        hi_coeff = 0.5 / (hi_af["total_exptime"] / sci_exp)
        combined_dark = linear_combination([lo_dark, hi_dark], [lo_coeff, hi_coeff])
# TO DO, plots?
# TO DO, what is excluded rows?
        excluded_rows = [2,3,4,7,8]
        for i in range(binned_sci.shape[0]):
            if i not in excluded_rows:
                plt.plot(binned_sci[i], label=i)
        plt.legend()
        plt.savefig("rows.png")
        plt.clf()
        val_c = c_stat(combined_dark, binned_sci, excluded_rows)

# TO DO, always hardcode this?
        x0 = [0.007, 0.005]
        res = minimize(fun_opt, x0, method="Nelder-Mead", tol=1e-6, 
                       args=([lo_dark, hi_dark], binned_sci, excluded_rows))
        combined_dark1 = linear_combination([lo_dark, hi_dark], res.x)
        plt.imshow(combined_dark1, aspect="auto", origin="lower")
        plt.savefig("combined_dark.png")
        plt.clf()

        sci_row = 9
        plt.plot(binned_sci[9], label="Outside of sci extr")
        plt.plot(combined_dark1[9], label="Predicted dark level")
        plt.xlabel("X")
        plt.ylabel("Number of photons")
        plt.legend()
        plt.savefig("predicted_dark.png")
        plt.clf()

    lo_af.close()
    hi_af.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datadir",
                        help="Path to science corrtags")
    parser.add_argument("--lo", dest="lo_darkname",
                        help="Name of low activity superdark")
    parser.add_argument("--hi", dest="hi_darkname",
                        help="Name of high activity superdark")
    args = parser.parse_args()
    corrtags = glob.glob(os.path.join(args.datadir, "*corrtag*fits"))

    main(corrtags, args.lo_darkname, args.hi_darkname)

