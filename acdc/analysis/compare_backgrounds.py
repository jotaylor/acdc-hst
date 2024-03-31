import os
import itertools
import numpy as np
import argparse
from matplotlib import pyplot as plt
dirname = os.path.dirname(__file__)
stylesheet = os.path.join(dirname, "niceplot.mplstyle")
plt.style.use(stylesheet)
import glob
from astropy.io import fits
from scipy import stats
import asdf

from acdc.superdark.cos_fuv_superdark import Superdark
from acdc.utils.utils import get_binning_pars, bin_coords

LOCAL_REFDIR = "/grp/hst/cdbs/lref"
LOCAL_DARKDIR = "/astro/sveash/cos_dark/final_superdarks"
PHA_INCLUSIVE = [2, 23]
PHA_INCL_EXCL = [2, 24]
# Colorblind-safe palette below
COLORS = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#a6cee3",
          "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
          "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#a6cee3"]


def sum_bkg_region(data, bkg_lo, bkg_hi):
    assert np.shape(data)[1] == len(bkg_lo) == len(bkg_hi), "Data X-dimension must match limits of background regions"
    datasum = []
    for i in range(len(bkg_lo)):
        colsum = np.sum(data[bkg_lo[i]:bkg_hi[i], i])
        datasum.append(colsum)
    datasum = np.array(datasum)
    
    return datasum

def smooth_array(data, smoothsize=100):
    bin_edges = np.arange(0, len(data), smoothsize)
    binned_y, edges, binnumber = stats.binned_statistic(np.arange(len(data)),
                                                      data,
                                                      statistic="mean",
                                                      bins=bin_edges)
    binned_x = bin_edges[:-1] + (smoothsize/2.)
    return binned_y, binned_x

def get_data(flt_custom, flt_default, superdark_custom, dark_dir=LOCAL_DARKDIR, ref_dir=LOCAL_REFDIR):
    data_custom = fits.getdata(flt_custom)
    data_default = fits.getdata(flt_default)
    segment = fits.getval(flt_custom, "segment")
    cenwave = fits.getval(flt_custom, "cenwave")
    hv = fits.getval(flt_custom, f"HVLEVEL{segment[-1]}", 1)
    xtractab_file0 = fits.getval(flt_custom, "xtractab")
    xtractab_file1 = xtractab_file0.split("$")[1]
    xtractab_file = os.path.join(ref_dir, xtractab_file1)

    xtract = fits.getdata(xtractab_file)
    inds = np.where((xtract["segment"] == segment) &
                    (xtract["cenwave"] == cenwave) & 
                    (xtract["aperture"] == "PSA"))
    xtract_data = xtract[inds]
    b1 = xtract_data["b_bkg1"]
    b2 = xtract_data["b_bkg2"]
    height = xtract_data["height"]
    m = xtract_data["slope"]
    smoothsize = xtract_data["bwidth"]
    x = np.arange(0, 16384)
    bkg1_mid = m*x + b1
    bkg1_y1 = bkg1_mid + height/2
    bkg1_y0 = bkg1_mid - height/2
    bkg2_mid = m*x + b2
    bkg2_y1 = bkg2_mid + height/2
    bkg2_y0 = bkg2_mid - height/2
    bkg1_y0 = np.round(bkg1_y0)
    bkg1_y0 = bkg1_y0.astype(int)
    bkg2_y0 = np.round(bkg2_y0)
    bkg2_y0 = bkg2_y0.astype(int)
    bkg1_y1 = np.round(bkg1_y1)
    bkg1_y1 = bkg1_y1.astype(int)
    bkg2_y1 = np.round(bkg2_y1)
    bkg2_y1 = bkg2_y1.astype(int)

    superdark_lo = glob.glob(os.path.join(dark_dir, f"*{segment}*{hv}*quiescent*asdf"))
    superdark_lo = [x for x in superdark_lo if "binned" not in x]
    superdark_hi = glob.glob(os.path.join(dark_dir, f"*{segment}*{hv}*active*asdf"))
    superdark_hi = [x for x in superdark_hi if "binned" not in x]
    assert len(superdark_lo) != 0, \
        f"No quiescent superdarks found for {segment} {hv}" 
    assert len(superdark_hi) != 0, \
        f"No active superdarks found for {segment} {hv}" 
    assert len(superdark_lo) == 1, \
        f"Too many quiescent superdarks found for {segment} {hv}" 
    assert len(superdark_hi) == 1, \
        f"Too many active superdarks found for {segment} {hv}"
    superdark_lo = superdark_lo[0]
    superdark_hi = superdark_hi[0]
    S_lo = Superdark.from_asdf(superdark_lo) 
    S_lo.bin_superdark(1, 1, pha_bins=[1,31]) 
    #S_lo.bin_superdark(1, 1, pha_bins=PHA_INCL_EXCL) 
    S_hi = Superdark.from_asdf(superdark_hi) 
    S_hi.bin_superdark(1, 1, pha_bins=[1,31])
    #S_hi.bin_superdark(1, 1, pha_bins=PHA_INCL_EXCL)
    lo_exptime = S_lo.total_exptime
    hi_exptime = S_hi.total_exptime
    S_lo_flt0 = S_lo.superdarks[0] / lo_exptime
    S_hi_flt0 = S_hi.superdarks[0] / hi_exptime
    assert np.shape(S_lo_flt0) == np.shape(S_hi_flt0), \
        "Quiescent and active superdark dimensions do not match"
    xstart = S_lo.bin_xstart
    xend = S_lo.bin_xend
    ystart = S_lo.bin_ystart
    yend = S_lo.bin_yend
    S_lo_flt = np.zeros(16777216).reshape(1024, 16384)
    S_hi_flt = np.zeros(16777216).reshape(1024, 16384)
    S_lo_flt[ystart:yend, xstart:xend] = S_lo_flt0
    S_hi_flt[ystart:yend, xstart:xend] = S_hi_flt0

    bkg1_sum_lo = sum_bkg_region(S_lo_flt, bkg1_y0, bkg1_y1)
    bkg1_lo,_ = smooth_array(bkg1_sum_lo)
    bkg2_sum_lo = sum_bkg_region(S_lo_flt, bkg2_y0, bkg2_y1)
    bkg2_lo,_ = smooth_array(bkg2_sum_lo)
    bkg1_sum_hi = sum_bkg_region(S_hi_flt, bkg1_y0, bkg1_y1)
    bkg1_hi,_ = smooth_array(bkg1_sum_hi)
    bkg2_sum_hi = sum_bkg_region(S_hi_flt, bkg2_y0, bkg2_y1)
    bkg2_hi,_ = smooth_array(bkg2_sum_hi)

    bkg1_sum_custom = sum_bkg_region(data_custom, bkg1_y0, bkg1_y1)
    bkg1_custom,_ = smooth_array(bkg1_sum_custom)
    bkg2_sum_custom = sum_bkg_region(data_custom, bkg2_y0, bkg2_y1)
    bkg2_custom,_ = smooth_array(bkg2_sum_custom)
    bkg1_sum_default = sum_bkg_region(data_default, bkg1_y0, bkg1_y1)
    bkg1_default,_ = smooth_array(bkg1_sum_default)
    bkg2_sum_default = sum_bkg_region(data_default, bkg2_y0, bkg2_y1)
    bkg2_default, binned_x = smooth_array(bkg2_sum_default)
    
    af = asdf.open(superdark_custom)
    pha_str = f"pha{PHA_INCL_EXCL[0]}-{PHA_INCL_EXCL[1]}"
    dark_im = af.tree[pha_str]
    dark_im_flt0 = dark_im / af.tree["total_exptime"]
    binning = get_binning_pars(af)
    dark_im_flt_perpixel = dark_im_flt0 / binning["bin_y"] / binning["bin_x"]
    dark_im_flt = np.zeros(16777216).reshape(1024, 16384)
    xs = np.arange(binning["xstart"], binning["xend"], binning["bin_x"])
    ys = np.arange(binning["ystart"], binning["yend"], binning["bin_y"])
    for i in range(len(xs)-1): 
        for j in range(len(ys)-1):
            dark_im_flt[ys[j]:ys[j+1], xs[i]:xs[i+1]] = dark_im_flt_perpixel[j,i]
    dark_x = np.arange(16384)
    
    bkg1_dark = sum_bkg_region(dark_im_flt, bkg1_y0, bkg1_y1)
    bkg2_dark = sum_bkg_region(dark_im_flt, bkg2_y0, bkg2_y1)
    
    return binned_x, bkg1_default, bkg2_default, bkg1_custom, bkg2_custom, \
        bkg1_lo, bkg2_lo, bkg1_hi, bkg2_hi, bkg1_dark, bkg2_dark, dark_x
    #return np.arange(16384), bkg1_sum_default, bkg2_sum_default, bkg1_sum_custom, bkg2_sum_custom, bkg1_sum_lo, bkg2_sum_lo, bkg1_sum_hi, bkg2_sum_hi


def plot_data(x, bkg1_default, bkg2_default, bkg1_custom, bkg2_custom, bkg1_lo, 
              bkg2_lo, bkg1_hi, bkg2_hi, bkg1_dark, bkg2_dark, dark_x, targname):
    fig, axes0 = plt.subplots(3, 1, figsize=(8.5, 11))
    axes = axes0.flatten()
    axes[0].hlines(y=0, xmin=-1000, xmax=17384, colors="lightgrey", linestyles="--")
    axes[1].hlines(y=0, xmin=-1000, xmax=17384, colors="lightgrey", linestyles="--")
    axes[2].hlines(y=0, xmin=-1000, xmax=17384, colors="lightgrey", linestyles="--")
    
    bkg1 = {"default": bkg1_default, "custom": bkg1_custom, "lo": bkg1_lo, "hi": bkg1_hi, "superdark": bkg1_dark}
    bkg2 = {"default": bkg2_default, "custom": bkg2_custom, "lo": bkg2_lo, "hi": bkg2_hi, "superdark": bkg2_dark}
    av = {"default": (bkg1_default + bkg2_default) / 2.,
          "custom": (bkg1_custom + bkg2_custom) / 2.,
          "lo": (bkg1_lo + bkg2_lo) / 2.,
          "hi": (bkg1_hi + bkg2_hi) / 2.,
          "superdark": (bkg1_dark + bkg2_dark) / 2.}
    dicts = [bkg1, bkg2, av]

    titles = ["Background 1 (Lower)", "Background 2 (Upper)", "Average"]
    for i in range(3):
        ax = axes[i]
        d = dicts[i]
        ax.set_ylabel("Counts/s")
        ax.plot(x, d["default"],   COLORS[0], alpha=0.8, label="Default FLT")
        ax.plot(x, d["custom"],    COLORS[1], alpha=0.8, label="Custom FLT")
        ax.plot(x, d["lo"],        COLORS[2], alpha=0.8, label="Lo Dark")
        ax.plot(x, d["hi"],        COLORS[3], alpha=0.8, label="Hi Dark")
        ax.plot(dark_x, d["superdark"], COLORS[4], alpha=0.8, label="Custom Dark")
        ax.set_ylim(-.0001, 0.00035)
        ax.set_xlim(-500, 16884)
        ax.set_title(titles[i])
        ax.legend()
    figname = f"{targname}_flt_comp.png"
    fig.savefig(figname)
    print(f"Wrote {figname}")


def compare_backgrounds(flt_default, flt_custom, targname, superdark_custom):
    binned_x, bkg1_default, bkg2_default, bkg1_custom, bkg2_custom, bkg1_lo, \
        bkg2_lo, bkg1_hi, bkg2_hi, bkg1_dark, bkg2_dark, dark_x = get_data(flt_custom, flt_default, superdark_custom)
    plot_data(binned_x, bkg1_default, bkg2_default, bkg1_custom, bkg2_custom, \
        bkg1_lo, bkg2_lo, bkg1_hi, bkg2_hi, bkg1_dark, bkg2_dark, dark_x, targname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--default",
                        help="Name of default CalCOS FLT file")
    parser.add_argument("-c", "--custom",
                        help="Name of custom-corrected FLT file")
    parser.add_argument("-t", "--targname",
                        help="Name of target")
    parser.add_argument("-s", "--superdark_custom",
                        help="Name of combined superdark from ACDC")
    args = parser.parse_args()
    compare_backgrounds(args.default, args.custom, args.targname, args.superdark_custom)

