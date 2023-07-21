import datetime
import numpy as np
import pandas as pd
from acdc.database.calculate_dark import get_aperture_region

def sql_to_df(sql_results, returncols):
    # Convert query results (list of attributes) to pandas dataframe
#    d = {col: [getattr(x, col) for x in results] for col in returncols}
    d = {}
    for col in returncols:
        try:
            d[col] = [getattr(x, col) for x in sql_results]
        except:
            pass
    df = pd.DataFrame(d)

    return df

def timefunc(func):  # IN USE
    """
    Decorator to wrap and time functions.

    Parameters
    ----------
    func: func
        Function to be timed

    Returns
    -------
    wrapper: func
        Timed version of the function
    """

    def wrapper(*args, **kw):
        t1 = datetime.datetime.now()
        result = func(*args, **kw)
        t2 = datetime.datetime.now()

        print(f"{func.__name__} executed in {t2 - t1}")

        return result

    return wrapper

def get_binning_pars(af):
    """
    For a given superdark, return the binning information in the spatial     
    directions and PHA.                                       
    """                                                       
    keys = ["bin_pha", "bin_x", "bin_y", "xstart", "xend", "ystart", "yend", 
            "phastart", "phaend"]                             
    binning = {}                                              
    for k in keys:                                            
        binning[k] = af[k]                                    
                                                              
    return binning                                            

def bin_coords(xs, ys, bin_x, bin_y, xstart=0, ystart=0, make_int=False, as_float=False):
    """
    Given a list of coordinates in X & Y, transform them into the superdark's
    binned (and possibly offset) coordinate system.
    """
    
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if not isinstance(ys, np.ndarray):
        ys = np.array(ys)
    if as_float is True:
        xsnew = (xs - xstart) / bin_x
        ysnew = (ys - ystart) / bin_y
    else:
        xsnew = (xs - xstart) // bin_x
        ysnew = (ys - ystart) // bin_y
    if make_int is True:
        xsnew = xsnew.astype(int)
        ysnew = ysnew.astype(int)

    return xsnew, ysnew

def unbin_coords(xs, ys, bin_x, bin_y, xstart=0, ystart=0):
    """
    Given a list of binned coordinates in X & Y, transform them into the the
    unbinned, native coordinate system.
    """
    
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if not isinstance(ys, np.ndarray):
        ys = np.array(ys)
    xsnew0 = (xs*bin_x) + xstart
    ysnew0 = (ys*bin_y) + ystart
    xsnew1 = xsnew0 + bin_x - 1
    ysnew1 = ysnew0 + bin_y - 1
    return (xsnew0, xsnew1), (ysnew0, ysnew1)

def unbin_image(binned_im, bin_x, bin_y, xstart=0, ystart=0, xend=16384, yend=1024):
    im_perpixel = binned_im / bin_x / bin_y
    unbinned_im = np.zeros(16777216).reshape(1024, 16384)
    xs = np.arange(xstart, xend, bin_x)
    ys = np.arange(ystart, yend, bin_y)
    for i in range(len(xs)-1):
        for j in range(len(ys)-1):
            unbinned_im[ys[j]:ys[j+1], xs[i]:xs[i+1]] = im_perpixel[j,i]
    return unbinned_im

def get_psa_wca(segment, cenwave, lp, binning, pad_psa=[0,0], pad_wca=[0,0]):
    """
    Determine the rows that correspond to the PSA and WCA apertures for a 
    given segment, cenwave, and lifetime position. Return the row indices in 
    the binned superdark coordinate system.

    Arguments:
        segment (str): Get PSA/WCA regions for this segment.
        cenwave (str): Get PSA/WCA regions for this cenwave.
        lp (str): Get PSA/WCA regions for this lifetime position.
        binning (dict): Dictionary that describes the binning information
            in both the spatial and PHA dimensions.
    Returns:
        excluded_rows (array): List of rows to exclude. Corresponds
            to the PSA and WCA regions.
        apertures (dict): The excluded rows for each aperture, PSA and WCA.
    """

    excluded_rows = np.array(())
    apertures = {"PSA": None, "WCA": None}
    pad = {"PSA": pad_psa, "WCA": pad_wca}
    for aperture in apertures:
        aperture_regions = get_aperture_region(cenwave=cenwave, aperture=aperture, segments=[segment], life_adj=[lp])
        box = aperture_regions[segment][f"lp{lp}_{aperture.lower()}_{cenwave}"]
        xmin0, xmax0, ymin0, ymax0 = box                            
        xsnew, ysnew = bin_coords(np.array([xmin0, xmax0]), np.array([ymin0, ymax0]), binning["bin_x"], binning["bin_y"], binning["xstart"], binning["ystart"], as_float=True, make_int=False)
#        ap_xmin, ap_xmax = xsnew # We don't need x coords
        ap_ymin, ap_ymax = ysnew
        ap_ymin = int(np.floor(ap_ymin))
        ap_ymax = int(np.ceil(ap_ymax))
        ap_ymin -= pad[aperture][0]
        ap_ymax += pad[aperture][1]
        rows = np.arange(ap_ymin, ap_ymax)
        apertures[aperture] = rows
        excluded_rows = np.concatenate((excluded_rows, rows))
    excluded_rows = excluded_rows.astype(int)
    return excluded_rows, apertures

