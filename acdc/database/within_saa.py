import numpy as np
import math

# Model 31 from https://github.com/spacetelescope/costools/blob/master/costools/saamodel.py
COS_FUV_MODEL = [
	(-28.3,         14.0),
	(-27.5,         15.0),
	(-26.1,         13.0),
	(-19.8,          1.5),
	( -9.6,        341.0),
	( -7.6,        330.4),
	( -6.0,        318.8),
	( -7.9,        297.2),
	(-12.0,        286.1),
	(-17.1,        279.9),
	(-20.3,        277.5),
	(-23.5,        276.5),
	(-26.0,        276.4),
	(-28.6,        276.7)]

DEGtoRAD = math.pi / 180.
TWOPI = 2. * math.pi
# This cutoff is based on current models in saamodel.py; this is used when
# finding the middle of the SAA region (middle_SAA).
SAA_LONGITUDE_CUTOFF = 200.



def testWithinSAA(hst, vertices, middle_SAA):
    """Test whether HST is within the polygon for an SAA contour.
    
    Args:
        hst (array_like): Unit vector pointing from the center of the Earth toward the
            location of HST at a particular time.
        vertices (array_like, shape (nvertices,3)): vertices[i] is a unit vector 
            from the center of the Earth toward
            vertex number i of a polygon that defines one of the SAA contour.
        middle_SAA (array_like): Unit vector from the center of the Earth toward a point near the
            middle of the SAA region.  This is for making a quick check that
            hst is close enough to the SAA contour to be worth making a
            detailed check.
    
    Returns:
        bool: True if hst is within the SAA contour defined by vertices.
    """

    # This test is primarily to exclude points that are diametrically
    # opposite to the SAA contour (because this would not be caught by
    # the code below!), but it should also save unnecessary arithmetic
    # most of the time.
    if np.dot(hst, middle_SAA) < 0.:
        return False

    nvertices = len(vertices)

    sin_lat_hst = hst[2]
    cos_lat_hst = math.sqrt(1. - sin_lat_hst**2)
    cos_long_hst = hst[0] / cos_lat_hst
    sin_long_hst = hst[1] / cos_lat_hst

    # vertices rotated to put hst in the x-z plane
    v_rot = vertices.copy()
    v_rot[:,0] =  vertices[:,0] * cos_long_hst + vertices[:,1] * sin_long_hst
    v_rot[:,1] = -vertices[:,0] * sin_long_hst + vertices[:,1] * cos_long_hst

    # v_rot rotated to put hst on the x axis
    v_rotrot = v_rot.copy()
    v_rotrot[:,0] =  v_rot[:,0] * cos_lat_hst + v_rot[:,2] * sin_lat_hst
    v_rotrot[:,2] = -v_rot[:,0] * sin_lat_hst + v_rot[:,2] * cos_lat_hst

    azimuth = np.arctan2(v_rotrot[:,2], v_rotrot[:,1])
    azimuth = np.where(azimuth < 0., azimuth + TWOPI, azimuth)

    delta_az = azimuth[1:] - azimuth[0:-1]
    delta_az = np.where(delta_az < -np.pi, delta_az + TWOPI, delta_az)
    delta_az = np.where(delta_az > np.pi, delta_az - TWOPI, delta_az)

    sum_delta_az = delta_az.sum()

    return not (sum_delta_az < 0.1 and sum_delta_az > -0.1)


def saaFilter(longitude_col, latitude_col, model=COS_FUV_MODEL):
    """Flag within the specified SAA contour as bad.
    
    Args:
        model (int): The SAA model number.  Currently these range from 2 to 32
            inclusive.  (Models 0 and 1 are radio frequence interference
            contours.)
    
    Returns:
        flag (array_like): This is a boolean array, one element for each row of the
            TIMELINE table.  True means that HST was within the SAA
            contour (specified by model) at the time corresponding to
            the TIMELINE row.
    """

    nelem = len(longitude_col)

    flag = np.zeros(nelem, dtype=np.bool8)

    model_vertices = model
    model_vertices.append(model_vertices[0])        # make a closed loop
    nvertices = len(model_vertices)
    # will be unit vectors from center of Earth pointing toward vertices
    vertices = np.zeros((nvertices, 3), dtype=np.float64)
    minmax_long = [720., -360.]     # will be minimum, maximum longitudes
    minmax_lat = [90., -90.]        # will be minimum, maximum latitudes
    for i in range(nvertices):
        (latitude, longitude) = model_vertices[i]
        vertices[i] = toRect(longitude, latitude)   # change the order
        if longitude < SAA_LONGITUDE_CUTOFF:
            longitude += 360.
        minmax_long[0] = min(longitude, minmax_long[0])
        minmax_long[1] = max(longitude, minmax_long[1])
        minmax_lat[0] = min(latitude, minmax_lat[0])
        minmax_lat[1] = max(latitude, minmax_lat[1])
    middle_long = (minmax_long[0] + minmax_long[1]) / 2.
    middle_lat = (minmax_lat[0] + minmax_lat[1]) / 2.
    middle_SAA = toRect(middle_long, middle_lat)

    # for each row in TIMELINE table
    for k in range(nelem):
        hst = toRect(longitude_col[k], latitude_col[k])
        flag[k] = testWithinSAA(hst, vertices, middle_SAA)

    return flag


def toRect(longitude, latitude):
    """Convert longitude and latitude to rectangular coordinates.
    
    Args:
        longitude (float): longitude in degrees.
        latitude (float): latitude in degrees.
    
    Returns:
        rect (array-like): Unit vector in rectangular coordinates.
    """

    rect = np.array([1.0, 0.0, 0.0], dtype=np.float64)
    longitude *= DEGtoRAD
    latitude *= DEGtoRAD

    rect[0] = math.cos(latitude) * math.cos(longitude)
    rect[1] = math.cos(latitude) * math.sin(longitude)
    rect[2] = math.sin(latitude)

    return rect

