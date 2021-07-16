import argparse
import copy
import asdf
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import glob
import os

from predict_dark_level1 import bin_science, get_binning_pars

def main(corrtags, pred_noise_file, fact=2, outdir="."):
# fact changes resolution, 1=highest resolution, 8=lowest resolution. to keep at same input res., use 1
    pred_noise_af = asdf.open(pred_noise_file)
    pred_noise = pred_noise_af["3-29"]
    b = get_binning_pars(pred_noise_af)
    for item in corrtags:
        binned, nevents = bin_science(item, b, fact)
        
        hdulist = fits.open(item, mode="update")
        data = hdulist[1].data
        
        xstart = b["xstart"]
        xend = b["xend"]
        ystart = b["ystart"]
        yend = b["yend"]
        phastart = b["phastart"]
        phaend = b["phaend"]
        bin_x = b["bin_x"]
        bin_y = b["bin_y"]

        wwf = open("tst.txt", "w")
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
                        outstr = f"{data['epsilon'][i]}\t{delta_eps}\t{delta_eps_ind}\n"
                        data["epsilon"][i] -= delta_eps_ind
                        wwf.write(outstr)
        wwf.close()

        plt.imshow(logic, aspect="auto", origin="lower")
        plt.savefig("logic.png")
        plt.clf()
        
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
        hdulist[1].data = data
        outfile = os.path.join(outdir, f"corrected_{os.path.basename(item)}")
        hdulist.writeto(outfile, overwrite=True)
        print(f"Wrote {outfile}")

    pred_noise_af.close() 

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--datadir",
                        help="Path to science corrtags")
    parser.add_argument("-n", "--noise",
                        help="Predicted noise level file")
    parser.add_argument("-o", "--outdir", default=".",
                        help="Output directory")
    args = parser.parse_args()
    corrtags = glob.glob(os.path.join(args.datadir, "*corrtag*fits"))

    main(corrtags, args.noise, outdir=args.outdir)

