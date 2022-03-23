import datetime
import argparse
import asdf
from astropy.io import fits
from astropy.table import Table
import pandas as pd
import numpy as np
from calcos import ccos

from query_cos_dark import files_by_mjd

TESTING = False

# Re-binning method from 
# https://stackoverflow.com/questions/14916545/numpy-rebinning-a-2d-array
# This is here for reference, original definitions from calculate_dark.py
#if segment == "FUVA":
#    location = {"inner": (1260, 15119, 375, 660), 
#                "bottom": (1060, 15250, 296, 375),
#                "top": (1060, 15250, 660, 734), 
#                "left": (1060, 1260, 296, 734),
#                "right": (15119, 15250, 296, 734)}
#elif segment == "FUVB":
#    location = {"inner": (1000, 14990, 405, 740), 
#                "bottom": (809, 15182, 360, 405),
#                "top": (809, 15182, 740, 785), 
#                "left": (809, 1000, 360, 785),
#                "right": (14990, 15182, 360, 785)}

def make_clean_superdark(hv, segment, start_mjd, ndays=300, dayint=100, 
                         gsagtab="41g2040ol_gsag.fits", bin_pha=1, bin_x=1,
                         bin_y=1):
    # Inner region only, this divides evenly by resel bin sizes
    if segment == "FUVA":
        x0 = 1264
        x1 = 15112 
        y0 = 376
        y1 = 660
    else:
        x0 = 1000
        x1 = 14984 
        y0 = 405
        y1 = 739 
    b_y0 = y0//bin_y
    b_y1 = y1//bin_y
    b_x0 = x0//bin_x
    b_x1 = x1//bin_x

    superdark = np.zeros(1024*16384).reshape((1024, 16384))
    extfound = False
    extnum = 1
    hvkey = f"HVLEVEL{segment[-1]}"
    while extfound is False:
        try:
            gsag_hv = fits.getval(gsagtab, hvkey, extnum)
            gsag_seg = fits.getval(gsagtab, "segment", extnum)
        except KeyError:
            extnum += 1
            continue
        if gsag_hv == hv and gsag_seg == segment:
            extfound = True
        else:
            extnum += 1
    gsag = Table.read(gsagtab, format="fits", hdu=extnum)
    df0 = gsag.to_pandas()
    df0 = df0.rename(columns={"LX": "X0", "LY": "Y0"})
    df0["X1"] = df0["X0"] + df0["DX"]
    df0["Y1"] = df0["Y0"] + df0["DY"]
    df0["X0"] = df0["X0"]//bin_x
    df0["X1"] = df0["X1"]//bin_x
    df0["Y0"] = df0["Y0"]//bin_y
    df0["Y1"] = df0["Y1"]//bin_y
    df0 = df0.astype("int32")

    notfilled = True
    pha_images = {}
    start = start_mjd
    total_days = 0
    total_exptime = 0
    total_files = 0
    pha_range = np.arange(1, 31, bin_pha)
    if pha_range[-1] != 30:
        pha_range = np.concatenate( (pha_range, np.array([30])) )
    
    runstart = datetime.datetime.now()
    print("Start time: {}".format(runstart))
    while notfilled is True:
        total_days += dayint
        end = start + dayint
        notfilled = False
        print(f"Using darks from {start}-{end}, {total_days}/{ndays} days...")
        df = df0.loc[(df0["DATE"] > start) & (df0["DATE"] < end)]
        dark_df = files_by_mjd(start, end, segment=segment, hv=hv)
        darks = dark_df["fileloc"].values
        total_exptime += sum([fits.getval(x, "exptime", 1) for x in darks])
        total_files += len(darks)
        for i in range(len(pha_range)-1):
            sum_image = bin_corrtag(darks, phastart=pha_range[i], phaend=pha_range[i+1])
            tmp = sum_image.reshape(1024 // bin_y, bin_y, 16384 // bin_x, bin_x)
            binned = tmp.sum(axis=3).sum(axis=1)
            binned_inner = binned[b_y0:b_y1, b_x0:b_x1]
            for j in range(len(df)): 
                binned_inner[df.iloc[j]["Y0"]:df.iloc[j]["Y1"], df.iloc[j]["X0"]:df.iloc[j]["X1"]] = 999
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
        if total_days >= ndays:
            print("Not every pixel had events at every PHA")
            notfilled = False
    pha_images["xstart"] = b_x0 * bin_x
    pha_images["ystart"] = b_y0 * bin_y
    pha_images["phastart"] = pha_range[0]
    pha_images["xend"] = b_x1 * bin_x
    pha_images["yend"] = b_y1 * bin_y
    pha_images["phaend"] = pha_range[-1]
    pha_images["bin_x"] = bin_x
    pha_images["bin_y"] = bin_y
    pha_images["bin_pha"] = bin_pha
    pha_images["mjdstart"] = start_mjd
    end_mjd = start_mjd+total_days
    pha_images["mjdend"] = end_mjd
    pha_images["segment"] = segment
    pha_images["hv"] = hv
    pha_images["total_exptime"] = total_exptime
    pha_images["total_files"] = total_files
    af = asdf.AsdfFile(pha_images)
    today = runstart.strftime("%d%b%y")
    outfile = f"superdark_{today}_{segment}_{hv}_{start_mjd}_{end_mjd}.asdf"
    af.write_to(outfile)
    print(f"Wrote {outfile}")
    runend = datetime.datetime.now()
    print("End time: {}".format(runend))
    print("Total time: {}".format(runend-runstart))


def bin_corrtag(corrtag_list, phastart, phaend, xtype='XCORR', ytype='YCORR', sdqflags=0):
    """Bin one or more corrtag event lists into a 2D image.

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--hv", type=int, help="HV setting of interest")
    parser.add_argument("--segment", help="Segment of interest")
    parser.add_argument("--mjdstart", type=int, help="MJD start date")
    parser.add_argument("--ndays", type=int, default=100, 
                        help="Number of days beyond MJD start to gather data")
    parser.add_argument("--phastep", default=1, type=int,
                        help="Size of PHA binning")
    args = parser.parse_args() 
    make_clean_superdark(args.hv, args.segment, args.mjdstart, args.ndays,
                         bin_pha=args.phastep) 
