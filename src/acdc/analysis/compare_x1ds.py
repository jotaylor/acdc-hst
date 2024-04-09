import os
import datetime
import argparse
from matplotlib import pyplot as plt
dirname = os.path.dirname(__file__)
stylesheet = os.path.join(dirname, "niceplot.mplstyle")
plt.style.use(stylesheet)
import numpy as np
from astropy.io import fits
import glob
from scipy import stats
import pandas as pd

# Eventually put Svea's plots here

def compare_x1ds(custom_coadd_txt, default_coadd_txt, targ=None, outfile=None, 
                 save=True):

    custom = pd.read_csv(custom_coadd_txt, names=["i", "wl", "flux", "err"], 
                            header=None, delim_whitespace=True)   
    default = pd.read_csv(default_coadd_txt, names=["i", "wl", "flux", "err"], 
                            header=None, delim_whitespace=True)
    
    fig, axes0 = plt.subplots(3, 1, figsize=(20, 16), sharex=True)
    axes = axes0.flatten()
    
    axes[0].plot(custom["wl"], custom["flux"], "crimson", alpha=0.8, label="Custom")
    axes[0].plot(default["wl"], default["flux"], "royalblue", alpha=0.8, label="Default")
    axes[0].set_xlabel("Wavelength [A]")
    axes[0].set_ylabel("Flux [ergs/s/cm^2/A]")
    axes[0].legend(loc="lower right")

    axes[1].plot(custom["wl"], custom["err"], "crimson", alpha=0.8, label="Custom")
    axes[1].plot(default["wl"], default["err"], "royalblue", alpha=0.8, label="Default")
    axes[1].set_xlabel("Wavelength [A]")
    axes[1].set_ylabel("Error [ergs/s/cm^2/A]")
    axes[1].legend(loc="lower right")

    axes[2].plot(custom["wl"], custom["flux"]/custom["err"], "crimson", 
                 alpha=0.8, label="Custom")
    axes[2].plot(default["wl"], default["flux"]/default["err"], "royalblue", 
                 alpha=0.8, label="Default")
    axes[2].set_xlabel("Wavelength [A]")
    axes[2].set_ylabel("S/N")
    axes[2].legend(loc="lower right")

    fig.suptitle(targ)
    plt.tight_layout()

    if save is True:
        if outfile is None:
            if targ is not None:
                outfile = f"{targ}_x1d_comp.png"
            else:
                now = datetime.datetime.now()
                outfile = f"{now.strftime('%d%b%Y_%H%M%S')}_x1d_comp.png"
        plt.savefig(outfile, )
        print(f"Saved {outfile}")
    else:
        print("\nClose figure when finished")
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--custom",
                        help="Coadd txt file for custom dark-corrected data")
    parser.add_argument("-d", "--default",
                        help="Coadd txt file for default dark-corrected data")
    parser.add_argument("-t", "--targ", default=None,
                        help="Name of target being compared")
    parser.add_argument("-o", "--outfile", default=None,
                        help="Name of output figure")
    parser.add_argument("--interact", default=False, action="store_true",
                        help="If True, open plot interactive window and do not save figure")
    args = parser.parse_args()
    save = not args.interact

    compare_x1ds(args.custom, args.default, args.targ, args.outfile, save)
