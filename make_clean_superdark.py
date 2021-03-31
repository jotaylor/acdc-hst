import asdf
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import numpy as np
from calcos import ccos

from query_darks import files_by_mjd

# Re-binning method from 
# https://stackoverflow.com/questions/14916545/numpy-rebinning-a-2d-array
TESTING = False
HV = 167
SEGMENT = "FUVA"
START_DEC = 2017.
if TESTING is True:
    START_MJD = 57100
else:
    START_MJD = 57754
BIN_X = 8
BIN_Y = 2

# This is here for reference, these are the original definitions from calculate_dark.py
#if segment == "FUVA":
#    location = {"inner": (1260, 15119, 375, 660), "bottom": (1060, 15250, 296, 375),
#                "top": (1060, 15250, 660, 734), "left": (1060, 1260, 296, 734),
#                "right": (15119, 15250, 296, 734)}
#elif segment == "FUVB":
#    location = {"inner": (1000, 14990, 405, 740), "bottom": (809, 15182, 360, 405),
#                "top": (809, 15182, 740, 785), "left": (809, 1000, 360, 785),
#                "right": (14990, 15182, 360, 785)}

def make_clean_superdark(hv=HV, segment=SEGMENT, start_mjd=START_MJD, 
                         ndays=100, gsagtab="41g2040ol_gsag.fits",
                         pha_step=1):
    # Inner region only, this divides evenly by BIN_X and BIN_Y
    if segment == "FUVA":
        if TESTING is True:
            x0 = 8000
            x1 = 8024 
            y0 = 436
            y1 = 442
        else:
            x0 = 1264
            x1 = 15112 
            y0 = 376
            y1 = 660
    else:
        x0 = 1000
        x1 = 14984 
        y0 = 405
        y1 = 739 
#    x = (x1-x0)/BIN_X
#    y = (y1-y0)/BIN_Y
#    superdark = np.zeros(x*y).reshape((y,x))

    superdark = np.zeros(1024*16384).reshape((1024, 16384))
    gsag = Table.read(gsagtab, format="fits", hdu=28)
    df0 = gsag.to_pandas()
    df0 = df0.rename(columns={"LX": "X0", "LY": "Y0"})
    df0["X1"] = df0["X0"] + df0["DX"]
    df0["Y1"] = df0["Y0"] + df0["DY"]
    df0["Y0_INVERT"] = 1024-df0["Y1"]
    df0["Y1_INVERT"] = 1024-df0["Y0"]
    df0["X0"] = df0["X0"]//8
    df0["X1"] = df0["X1"]//8
    df0["Y0"] = df0["Y0"]//2
    df0["Y1"] = df0["Y1"]//2
    df0 = df0.astype("int32")

    notfilled = True
    pha_images = {}
    start = start_mjd
    total_days = 0
    total_exptime = 0
    pha_range = np.arange(1, 31, pha_step)
    if pha_range[-1] != 30:
        pha_range = np.concatenate( (pha_range, np.array([30])) )
    
    while notfilled is True:
        total_days += ndays
        end = start + ndays
        notfilled = False
        print(total_days)
        df = df0.loc[(df0["DATE"] > start) & (df0["DATE"] < end)]
        dark_df = files_by_mjd(start, end, segment=segment, hv=hv)
        darks = dark_df["fileloc"].values
        total_exptime += sum([fits.getval(x, "exptime", 1) for x in darks])
        for i in range(len(pha_range)-1):
            sum_image = bin_corrtag(darks, phastart=pha_range[i], phaend=pha_range[i+1])
#            inner_image = sum_image[(1024-y1):(1024-y0), x0:x1]
            tmp = sum_image.reshape(1024 // 2, 2, 16384 // 8, 8)
            binned = tmp.sum(axis=3).sum(axis=1)
            b_y0 = (1024-y1)//BIN_Y
            b_y1 = (1024-y0)//BIN_Y
            b_x0 = (16384-x1)//BIN_X
            b_x1 = (16384-x0)//BIN_X
            binned_inner = binned[b_y0:b_y1, b_x0:b_x1]
            for j in range(len(df)): 
                binned_inner[df.iloc[j]["Y0_INVERT"]:df.iloc[j]["Y1_INVERT"], df.iloc[j]["X0"]:df.iloc[j]["X1"]] = 999
            zeros = np.where(binned_inner == 0)
            if len(zeros[0]) != 0:
                notfilled = True
                if TESTING is True:
                    break
            key = f"pha{pha_range[i]}-{pha_range[i+1]}"
            if i in pha_images:
                pha_images[key] += binned_inner
            else:
                pha_images[key] = binned_inner
        start = end
        if total_days >= 300:
            notfilled = False
    pha_images["xstart"] = x0
    pha_images["ystart"] = y0
    pha_images["xend"] = x1
    pha_images["yend"] = y1
    pha_images["mjdstart"] = start_mjd
    pha_images["mjdend"] = start_mjd+total_days
    pha_images["segment"] = segment
    pha_images["hv"] = hv
    pha_images["total_exptime"] = total_exptime
    af = asdf.AsdfFile(pha_images)
    outfile = f"superdark_{segment}_{hv}.asdf"
    af.write_to(outfile)
    print(f"Wrote {outfile}")


def bin_corrtag(corrtag_list, phastart, phaend, xtype='XCORR', ytype='YCORR', sdqflags=0):
    """Bin corrtag event lists into a 2D image.

    Modifed from /user/jotaylor/git/jotools/utils.py

    Parameters
    ----------
    corrtag_list : str
        list of datasets to combine into an image
    xtype : str, optional
        X-coordinate type
    ytype : str, optional
        Y-coordinate type

    Returns
    -------
    image : np.ndarray
        2D image of the corrtag
    """

    final_image = np.zeros((1024, 16384)).astype(np.float32)

    for filename in corrtag_list:
        image = np.zeros((1024, 16384)).astype(np.float32)
        events = fits.getdata(filename, 1)
        inds = np.where((events["pha"] >= phastart) & (events["pha"] <= phaend))
        events = events[inds]
        if not len(inds[0]):
            print(f"No counts for PHA={pha}")
            return image
        
        # Call for this is x_values, y_values, image to bin to, offset in x
        # ccos.binevents(x, y, array, x_offset, dq, sdqflags, epsilon)
        ccos.binevents(events[xtype].astype(np.float32),
                       events[ytype].astype(np.float32),
                       image,
                       0,
                       events['dq'],
                       sdqflags)
        
        final_image += image
    
    return final_image


if __name__ == "__main__":
    make_clean_superdark() 
