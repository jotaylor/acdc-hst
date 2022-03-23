import argparse
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import fits
import glob
from scipy import stats
import pandas as pd

# Eventually put Svea's plots here

def compare_x1ds(targ, save=True):
#def compare_x1ds(custom_txt, default_txt, targ, save=True):
    custom_txt = glob.glob(f"custom_012722_{targ}*txt")[0]
    #custom_txt = glob.glob(f"custom_011222_{targ}*txt")[0]
    default_txt = glob.glob(f"default_011222_{targ}*txt")[0]
    custom = pd.read_csv(custom_txt, names=["i", "wl", "flux", "err"], 
                            header=None, delim_whitespace=True)   
    default = pd.read_csv(default_txt, names=["i", "wl", "flux", "err"], 
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
        outfile = f"{targ}_x1d_comp.png"
        plt.savefig(outfile, )
        print(f"Saved {outfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
#    parser.add_argument(dest="custom_txt")
#    parser.add_argument(dest="default_txt")
    parser.add_argument(dest="targ")
    args = parser.parse_args()

#    compare_x1ds(args.custom_txt, args.default_txt, args.targ)
    compare_x1ds(args.targ)
