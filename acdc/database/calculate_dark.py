import datetime
import os
import glob
import numpy as np
import pandas as pd
from ftplib import FTP
from astropy.io import fits
from astropy.time import Time
from astropy.io import ascii

__author__ = "Jo Taylor"
__email__ = "jotaylor@stsci.edu"

def measure_darkrate(filename, apertures=None):
    """
    For an input dark dataset, record exposure information including
    observation time, observatory latitude & longitude. Measure dark rate at
    specified regions of the COS FUV detector at a given time interval in order
    to accumulate enough counts for a measurement.
    (Taken from the STScI COS team"s dark monitor, dark_monitor.py, and modified. 
    Original authors include Justin Ely, Mees Fix, Dzhuliya Dashtamirova.)
    
    Args:
        filename (str): Name of dark dataset.
        apertures (dict): Output from `get_aperture_region()`. This is passed in to avoid
            running `get_aperture_region()` many times within this function.
   
    Returns:
        dark_df (:obj:`pandas.dataframe`): Pandas dataframe with information for each
            dark measurement. 
    """

    hdulist = fits.open(filename)
    
    if apertures is None:
        apertures = get_aperture_region()

    timeline = hdulist["timeline"].data
    events = hdulist["events"].data
    segment = hdulist[0].header["segment"]
    rootname = hdulist[0].header["rootname"]

    if segment == "FUVA":
        location = {"inner": (1260, 15119, 375, 660), "bottom": (1060, 15250, 296, 375),
                    "top": (1060, 15250, 660, 734), "left": (1060, 1260, 296, 734),
                    "right": (15119, 15250, 296, 734)}
        location.update(apertures[segment])
    elif segment == "FUVB":
        location = {"inner": (1000, 14990, 405, 740), "bottom": (809, 15182, 360, 405),
                    "top": (809, 15182, 740, 785), "left": (809, 1000, 360, 785),
                    "right": (14990, 15182, 360, 785)}
        location.update(apertures[segment])
    pha = (2, 23)
    timestep = 25 # in units of seconds
    times = timeline["time"][::timestep].copy()
    events = hdulist["events"].data

    lat = timeline["latitude"][::timestep][:-1].copy()
    lon = timeline["longitude"][::timestep][:-1].copy()
    second_per_mjd = 1.15741e-5
    mjd_per_step = hdulist[1].header["EXPSTART"] + times.copy().astype(np.float64) * second_per_mjd
    mjd = mjd_per_step[:-1]                                                   
    t = Time(mjd, format="mjd") 
    mjdtime = [x.value for x in t]
    
    # Creates different tables, one for each region of detector
    info_dict = {"segment": segment, "rootname": rootname,"latitude": lat, 
                 "longitude": lon, "time": mjdtime, "pha": np.arange(32), 
                 "exptime": timestep}
    dark_dict = {}
    for x in location:
        region_area = (location.get(x)[1] - location.get(x)[0]) * (location.get(x)[3] - location.get(x)[2])
        # In reality this is not an issue since YCORR=YFULL for darks.
        if "location" == "apertures":
            coord_x = "XCORR"
            coord_y = "YFULL"
        else:
            coord_x = "XCORR"
            coord_y = "YCORR"
        counts_pha = []
        for phabin in range(32):
            index = np.where((events["PHA"] == phabin) &
                             (events[coord_x] > location.get(x)[0]) &
                             (events[coord_x] < location.get(x)[1]) &
                             (events[coord_y] > location.get(x)[2]) &
                             (events[coord_y] < location.get(x)[3]))
            counts = np.histogram(events[index]["time"], bins=times)[0]
            counts_pha.append(counts)
        counts_time = np.dstack(counts_pha)[0]
        region_d = {"darks": counts_time, "xcorr_min": location[x][0], 
                    "xcorr_max": location[x][1], "ycorr_min": location[x][2], 
                    "ycorr_max": location[x][3], "region_area": region_area}
        dark_dict[x] = region_d
    hdulist.close() 
    
    return info_dict, dark_dict

def get_solar_data(solardir):
    """
    Pull solar data files from NOAA website. Solar data is FTPd from NOAA and 
    written to text files.
    (Taken from the STScI COS team"s dark monitor, dark_monitor.py, and modified. 
    Original authors include Justin Ely, Mees Fix, Dzhuliya Dashtamirova.)

    Args:
        solardir (str): Directory to write the solar data files to.
    """

    ftp = FTP("ftp.swpc.noaa.gov")
    ftp.login()
    ftp.cwd("/pub/indices/old_indices/")

    for item in sorted(ftp.nlst()):
        if item.endswith("_DSD.txt"):
            year = int(item[:4])
            if year >= 2008:
                print("Retrieving: {}".format(item))
                destination = os.path.join(solardir, item)
                ftp.retrbinary("RETR {}".format(item),
                               open(destination, "wb").write)
                os.chmod(destination, 0o777)


def to_mjd(row):
    fulltime = row["time-tag"]+"-01 00:00:00"
    mjd = Time(fulltime, scale="utc", format="iso").mjd
    return mjd


def parse_solar_json(infile):
    """
    Parse new NOAA solar data. Taken from:
    https://www.swpc.noaa.gov/products/solar-cycle-progression
    which leads to:
    https://services.swpc.noaa.gov/json/solar-cycle/observed-solar-cycle-indices.json
    """
    df = pd.read_json(infile)
    df = df.drop(["ssn", "smoothed_ssn", "observed_swpc_ssn", "smoothed_swpc_ssn"], axis=1)
    df = df.loc[df["time-tag"] >= '2009-01']
    df["mjd"] = df.apply(to_mjd, axis=1) 
    df = df.drop(["time-tag"], axis=1)
    return df["mjd"].values, df["f10.7"].values


def parse_solar_files(files):
    """
    Parse solar data text files and return date and flux.
    (Taken from the STScI COS team"s dark monitor, dark_monitor.py, and modified. 
    Original authors include Justin Ely, Mees Fix, Dzhuliya Dashtamirova.)
    
    Args:
        files (str or array-like): If a string, a check is done to see if string 
            is a file or a directory to be globbed. If not a string, array-like
            is assumed.
    
    Returns:
        date (:obj:`numpy.ndarray`): MJD of each solar flux measurement.
        flux (:obj:`numpy.ndarray`): Solar flux measurements.
    """

    date = []
    flux = []
    if isinstance(files, str):
        if os.path.isfile(files):
            input_list = [files]
        else:
            input_list = glob.glob(os.path.join(files, "*DSD.txt"))
    else:
        input_list = files
    input_list.sort()

    for item in input_list:
#        print("Reading {}".format(item))

        # clean up Q4 files when year-long file exists
        if ("Q4_" in item) and os.path.exists(item.replace("Q4_", "_")):
            print("Removing duplicate observations: {}".format(item))
            os.remove(item)
            continue

        # astropy.ascii no longer returns an empty table for empty files
        # Throws IndexError, we will go around it if empty.
        try:
            data = ascii.read(item, data_start=1, comment="[#,:]")
        except IndexError:                                    
            continue                                          
                                                              
        for line in data:                                     
            line_date = Time("{}-{}-{} 00:00:00".format(line["col1"],       
                                                        line["col2"],       
                                                        line["col3"]),      
                             scale="utc", format="iso").mjd   
            line_flux = line[3]                               
                                                              
            if line_flux > 0:                                 
                date.append(line_date)                        
                flux.append(line_flux)                        
                                                              
    return np.array(date), np.array(flux)                     

def get_aperture_region(cenwave=1291, aperture="PSA", segments=["FUVA", "FUVB"], life_adj=[1, 2, 3, 4, 5, 6]):
    """
    Determine the extraction box for a given cenwave (1291 by default) and
    aperture (PSA by default).
    
    Args:
        cenwave (int): Cenwave to look up.
        aperture (str): Aperture to look up.

    Returns:
        cenwave (int): Cenwave setting.
        apertures (dict): Dictionary where each key is segment, and the value
            is another dictionary where each key is lifetime position,
            and the value is a tuple with (xmin, xmax, ymin, ymax) representing
            the cenwave extraction box for that segment and LP combo.
    """    

    if isinstance(life_adj, int):
        life_adj = [life_adj] 
    today0 = datetime.datetime.now()
    today = today0.strftime("%Y-%m-%d")

    if "CRDS_PATH" not in os.environ:
        if os.path.exists("/grp/crds/cache"):
            os.environ["CRDS_PATH"] = "/grp/crds/cache"
        else:
            raise AssertionError("CRDS_PATH environment variable must first be defined")
    os.environ["CRDS_SERVER_URL"] = "https://hst-crds.stsci.edu"
    import crds
    current_pmap = crds.get_default_context() 

    # The active area limits are taken from the COS BRFTAB x1u1459il_brf.fits
    # This file will almost certainly not be updated, so hardcoding is okay.
    aa_xcorr = {"FUVA": [1060, 15250], "FUVB": [809, 15182]}
    # aa_ycorr = {"FUVA": [296, 734], "FUVB": [360, 785]}

    apertures = {"FUVA": {}, "FUVB": {}}
    # For each LP, determine the appropriate xtractab as returned by CRDS on the fly.
    for segment in segments:
#        for life_adj, date_obs in zip([1, 2, 3, 4], ["2010-01-01", "2014-01-01", "2016-01-01", "2018-01-01"]):
        for lp in life_adj:
            crds_1dx = crds.getrecommendations(parameters={"INSTRUME": "COS", 
                                "DETECTOR": "FUV", "LIFE_ADJ": lp, 
                                "OBSTYPE": "SPECTROSCOPIC", 
                                "DATE-OBS": today, "TIME-OBS": "00:00:00",
                                "CENWAVE": cenwave},
                            reftypes=["xtractab"], context=current_pmap, observatory="hst")
            lp_1dx = os.path.join(os.environ["CRDS_PATH"], "references/hst/", 
                                  crds_1dx["xtractab"])
    
            data_1dx = fits.getdata(lp_1dx, 1)
            ind = np.where((data_1dx["cenwave"] == cenwave) & 
                           (data_1dx["segment"] == segment) & 
                           (data_1dx["aperture"] == aperture))
            aperture_data = data_1dx[ind]
    
            x = np.arange(16384)
            y_center = aperture_data["slope"][0] * x + aperture_data["b_spec"][0]
            y_upper = round(max(y_center + aperture_data["height"][0] / 2))
            y_lower = round(min(y_center - aperture_data["height"][0] / 2))
    
            lpkey = f"lp{lp}_{aperture.lower()}_{cenwave}"

            # This  matches the format already defined in measure_darkrate()
            # i.e. xmin, xmax, ymin, ymax 
            # VERY IMPORTANT: Y coords are in YFULL, X in XCORR
            apertures[segment][lpkey] = (aa_xcorr[segment][0], aa_xcorr[segment][1],
                                        int(y_lower), int(y_upper))
    return apertures
