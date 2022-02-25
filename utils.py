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

