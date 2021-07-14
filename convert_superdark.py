from matplotlib.backends.backend_pdf import PdfPages
import copy
import asdf
import numpy as np
import matplotlib.pyplot as plt
try:
    plt.style.use("niceplot.mplstyle")
except OSError:
    pass

def bin_superdark(superdark, bin_x=8, bin_y=2, bin_pha=26, 
                  phastart=3, phaend=28, outdir="."):

    with asdf.open(superdark) as af:
        # First bin by PHA
        binned = copy.deepcopy(af.tree[f"pha{phastart}-{phastart+1}"])
        binned_ims = {}
        counter = 1
        for i in range(phastart+1, phaend+2):
            if counter == bin_pha:
                counter = 0
                binned_ims[f"{i-bin_pha}-{i}"] = binned
                binned = copy.deepcopy(af.tree[f"pha{i}-{i+1}"])
            else:
                binned += af.tree[f"pha{i}-{i+1}"]
            counter += 1

        sh = np.shape(binned)
        xdim = sh[1]
        ydim = sh[0]
        b_x0 = 0
        b_x1 = (xdim // bin_x) * bin_x
        b_y0 = 0
        b_y1 = (ydim // bin_y) * bin_y

        pdffile = os.path.join(outdir, superdark.replace("asdf", "pdf"))
        pdf = PdfPages(pdffile)
        for k,binned in binned_ims.items():
            binned = binned[b_y0:b_y1, b_x0:b_x1]
            binned_sh = np.shape(binned)
            binned_xdim = binned_sh[1]
            binned_ydim = binned_sh[0]
            tmp = binned.reshape(binned_ydim // bin_y, bin_y, binned_xdim // bin_x, bin_x)
            binned = tmp.sum(axis=3).sum(axis=1)
            binned_ims[k] = binned
            spl = k.split("-")
            print(f"For PHAs {spl[0]} through {spl[1]}:")
            print(f"\tTotal number of events: {np.sum(binned)}")
            print(f"\tMinimum number of events in a binned pixel: {np.min(binned)}")
            print(f"\tMean number of events per binned pixel: {np.mean(binned):.1f}")
            print(f"\t  Standard deviation: {np.std(binned):.1f}")
            fig, ax = plt.subplots(figsize=(20,5))
            im = ax.imshow(binned, origin="lower", cmap="inferno", 
                      vmax=.8*np.max(binned))
            fig.colorbar(im)
            
            ax.set_title(f"{af['segment']} HV={af['hv']} MJD {af['mjdstart']}-{af['mjdend']} PHA {spl[0]}-{spl[1]} (binned inner region)")
            pdf.savefig(fig)
        pdf.close()
        print(f"Wrote {pdffile}")

        for k in af.keys():
            if k.startswith("pha"):
                continue
            binned_ims[k] = af[k]
        binned_ims["bin_x"] *= bin_x
        binned_ims["bin_y"] *= bin_y
        binned_ims["bin_pha"] *= bin_pha
        binned_ims["xend"] = b_x1 * bin_x + af["xstart"]
        binned_ims["yend"] = b_y1 * bin_y + af["ystart"]
        binned_ims["phastart"] = phastart
        binned_ims["phaend"] = phaend

    af = asdf.AsdfFile(binned_ims)
    outfile = f"binned_{superdark}"
    af.write_to(outfile)
    print(f"Wrote {outfile}")

    return binned_ims

