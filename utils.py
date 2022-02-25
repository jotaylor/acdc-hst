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

def bin_coords(xs, ys, bin_x, bin_y, xstart=0, ystart=0):
    """
    Given a list of coordinates in X & Y, transform them into the superdark's
    binned (and possibly offset) coordinate system.
    """
    xsnew = (xs - xstart) // bin_x
    ysnew = (ys - ystart) // bin_y
    return xsnew, ysnew