import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time

def measure_darkrate(filename=None):
    """
    Taken from the STScI COS team"s dark monitor, dark_monitor.py, and modified. 
    Original authors include Justin Ely, Mees Fix, Dzhuliya Dashtamirova.
    """

    hdulist = fits.open(filename)

    timeline = hdulist["timeline"].data
    events = hdulist["events"].data
    segment = hdulist[0].header["segment"]
    rootname = hdulist[0].header["rootname"]

    if segment == "FUVA": 
        location = {"inner": (1260, 15119, 375, 660), "bottom": (1060, 15250, 296, 375),
                    "top": (1060, 15250, 660, 734), "left": (1060, 1260, 296, 734),
                    "right": (15119, 15250, 296, 734)}
    elif segment == "FUVB":
        location = {"inner": (1000, 14990, 405, 740), "bottom": (809, 15182, 360, 405),
                    "top": (809, 15182, 740, 785), "left": (809, 1000, 360, 785),
                    "right": (14990, 15182, 360, 785)}
    pha = (2, 23)
    timestep = 25 # in units of seconds
    times = timeline["time"][::timestep].copy()
    events = hdulist["events"].data

    # Creates 5 different tables, one for each region of detector
    d = {}
    for x in location:
        region_area = (location.get(x)[1] - location.get(x)[0]) * (location.get(x)[3] - location.get(x)[2])
        index = np.where((events["PHA"] > pha[0]) &
                         (events["PHA"] < pha[1]) &
                         (events["XCORR"] > location.get(x)[0]) &
                         (events["XCORR"] < location.get(x)[1]) &
                         (events["YCORR"] > location.get(x)[2]) &
                         (events["YCORR"] < location.get(x)[3]))
        unfiltered_pha = np.where((events["XCORR"] > location.get(x)[0]) &
                                  (events["XCORR"] < location.get(x)[1]) &
                                  (events["YCORR"] > location.get(x)[2]) &
                                  (events["YCORR"] < location.get(x)[3]))
        counts = np.histogram(events[index]["time"], bins=times)[0]
        counts_unfiltered_pha = np.histogram(events[unfiltered_pha]["time"], bins=times)[0]
        lat = timeline["latitude"][::timestep][:-1].copy()
        lon = timeline["longitude"][::timestep][:-1].copy()
        second_per_mjd = 1.15741e-5
        mjd_per_step = hdulist[1].header["EXPSTART"] + times.copy().astype(np.float64) * second_per_mjd
        mjd = mjd_per_step[:-1]                                                   
        t = Time(mjd, format="mjd")                                               
        decyear = t.decimalyear                                                   
        counts = counts / region_area / timestep                                             
        counts_unfiltered_pha = counts_unfiltered_pha / region_area / timestep               
        d[x] = pd.DataFrame({"region": x, "segment": segment, "darkrate": counts, "time": decyear,
                            "xcorr_min": location.get(x)[0], "xcorr_max": location.get(x)[1],    
                             "ycorr_min": location.get(x)[2], "ycorr_max": location.get(x)[3],   
                             "longitude": lon, "latitude": lat, "rootname": rootname, 
                             "unfiltered_pha_counts": counts_unfiltered_pha})    
    # Combine all 5 tables for different regions                                  
    dark_df = pd.concat([d["inner"], d["bottom"], d["top"], d["left"], d["right"]])       
    # Flag data taken near SAA where flag=1 denotes data outside the SAA.
    dark_df["saa_flag"] = np.where(dark_df.eval("latitude > 10 or longitude < 260"), 1, 0)
    hdulist.close() 
    
    return dark_df

def dark_edges(df):
    """
    Takes a data frame of all edges and splits it up into a dict of separate
    dataframes, one for each segment and edge of the detector
    :param df: pandas dataframe of darks for all programs
    :return: a dictionary of dataframes split up by the segment and region
    """
    dark_edges_dict = {}
    seg = ["FUVA", "FUVB"]
    loc = ["inner", "bottom", "top", "left", "right"]

    for x in seg:
        for y in loc:
            dark_edges_dict[x + "_" + y] = (df[(df["segment"] == x) & (df["region"] == y)])
    return dark_edges_dict

