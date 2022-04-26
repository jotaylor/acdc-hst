import numpy as np
import pandas as pd

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

def bin_coords(xs, ys, bin_x, bin_y, xstart=0, ystart=0, make_int=False):
    """
    Given a list of coordinates in X & Y, transform them into the superdark's
    binned (and possibly offset) coordinate system.
    """
    
    if not isinstance(xs, np.ndarray):
        xs = np.array(xs)
    if not isinstance(ys, np.ndarray):
        ys = np.array(ys)
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

