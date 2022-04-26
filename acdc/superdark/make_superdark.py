import numpy as np
import argparse
from astropy.io import fits
import os
import glob

from acdc.database.query_cos_dark import files_by_mjd

def sum_data(mjdstart, mjdend):
    df = files_by_mjd(mjdstart, mjdend)

    total_time_a = 0
    total_time_b = 0
    sumdark_a = np.zeros(16384*1024).reshape(1024, 16384)
    sumdark_b = np.zeros(16384*1024).reshape(1024, 16384)
    for i in range(len(df)):
        fileloc = df.iloc[i]["fileloc"]
        segment = df.iloc[i]["segment"]
        filename = os.path.basename(fileloc)
        pid = os.path.dirname(fileloc).split("/")[-1]
        cosmofile = glob.glob(os.path.join("/grp/hst/cos2/cosmo/", pid, filename.replace("corrtag", "counts")+"*"))[0]
        print(cosmofile)
        exptime = fits.getval(cosmofile, "exptime", 1)
        data = fits.getdata(cosmofile)
        if segment == "FUVA":
            total_time_a += exptime
            sumdark_a += data
        else:
            total_time_b += exptime
            sumdark_b += data

    superdark_a = sumdark_a / total_time_a
    superdark_b = sumdark_b / total_time_b
    return superdark_a, superdark_b, total_time_a, total_time_b

def make_superdark(mjdstart, mjdend):
    superdark_a, superdark_b, total_time_a, total_time_b = sum_data(mjdstart, mjdend)
    
    hdr = fits.Header()
    hdr["MJDSTART"] = mjdstart
    hdr["MJDEND"] = mjdend
    primary = fits.PrimaryHDU(header=hdr)
    
    hdr1 = fits.Header()
    hdr1["EXTNAME"] = "SUPERDARK_FUVA"
    hdr1["EXPTIME"] = total_time_a
    sci1 = fits.ImageHDU(superdark_a, header=hdr1, name="SUPERDARK_FUVA")
    
    hdr2 = fits.Header()
    hdr2["EXTNAME"] = "SUPERDARK_FUVA"
    hdr2["EXPTIME"] = total_time_b
    sci2 = fits.ImageHDU(superdark_b, header=hdr2, name="SUPERDARK_FUVB")

    hdu = fits.HDUList([primary, sci1, sci2])
    outname = f"superdark_{mjdstart}_{mjdend}.fits"
    hdu.writeto(outname, overwrite=True)
    print(f"Wrote file {outname}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(dest="mjdstart", 
                        help="Starting MJD time to select files")
    parser.add_argument(dest="mjdend", 
                        help="Ending MJD time to select files")
    args = parser.parse_args()

    make_superdark(args.mjdstart, args.mjdend)
