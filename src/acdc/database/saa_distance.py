""" 
Taken from:
# https://stackoverflow.com/questions/43440813/shortest-great-circle-distance-between-a-point-and-a-polygon-on-a-sphere-globe
"""

from shapely.geometry import Polygon, mapping
import shapely
import numpy as np
import math
from numpy import linalg as LA

def Pairwise(iterable):
    """Iterate through an itertable returning adjacent pairs.
    
    Args:
        iterable (array-like): An iterable
    
    Returns:
        tuple: Pairs of sequential, adjacent entries from the iterable
    """
    it = iter(iterable)
    a = next(it, None)
    for b in it:
        yield (a, b)
        a = b


def LatLonToXYZ(lat,lon,radius):
    """Convert a latitude-longitude pair to 3D-Cartesian coordinates.

    Args:
        lat (float): Latitude in degrees
        lon (float): Longitude in degrees
        radius (float): Radius in arbitrary units
    
    Returns:
        tuple: x,y,z coordinates in arbitrary units
    """
    lat = np.radians(lat)
    lon = np.radians(lon)
    x = radius * np.cos(lon) * np.cos(lat)
    y = radius * np.sin(lon) * np.cos(lat)
    z = radius * np.sin(lat)
    return x,y,z


def XYZtoLatLon(x,y,z):
    """Convert 3D-Cartesian coordinates to a latitude-longitude pair.

    Args:
        x (float): x-coordinate in arbitrary units
        y (float): y-coordinate in arbitrary units
        z (float): z-coordinate in arbitrary units
    
    Returns:
        tuple: lat,lon pair in degrees
    """
    radius = np.sqrt(x*x+y*y+z*z)
    lat = np.degrees(np.arcsin(z/radius))
    lon = np.degrees(np.arctan2(y, x))   
    return lat,lon


def Haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the Great Circle distance on Earth between two latitude-longitude
    points
    Args:
        lat1 (float): Latitude of Point 1 in degrees
        lon1 (float): Longtiude of Point 1 in degrees
        lat2 (float): Latitude of Point 2 in degrees
        lon2 (float): Longtiude of Point 2 in degrees
    Returns:
        float: Distance between the two points in kilometres
    """
    Rearth = 6371
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    #Haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a)) 
    return Rearth*c


#http://stackoverflow.com/a/1302268/752843
def NearestPointOnGC(alat1,alon1,alat2,alon2,plat,plon):
    """
    Calculate the location of the nearest point on a Great Circle to a query point
    Args:
        lat1 (float): Latitude of start of arc in degrees
        lon1 (float): Longtiude of start of arc in degrees
        lat2 (float): Latitude of end of arc in degrees
        lon2 (float): Longtiude of end of arc in degrees
        plat (float): Latitude of query point in degrees
        plon (float): Longitude of query point in degrees
    Returns:
        tuple: A (lat,lon) pair in degrees of the closest point
    """
    Rearth    = 6371 #km
    #Convert everything to Cartesian coordinates
    a1 = np.array(LatLonToXYZ(alat1,alon1,Rearth))
    a2 = np.array(LatLonToXYZ(alat2,alon2,Rearth))
    p = np.array(LatLonToXYZ(plat, plon, Rearth))
    G = np.cross(a1,a2) #Plane of the Great Circle containing A and B
    F = np.cross(p,G)   #Plane perpendicular to G that passes through query pt
    T = np.cross(G,F)   #Vector marking the intersection of these planes
    T = Rearth*T/LA.norm(T) #Normalize to lie on the Great Circle
    tlat,tlon = XYZtoLatLon(*T)
    return tlat,tlon


def DistanceToGCArc(alat,alon,blat,blon,plat,plon):
    """
    Calculate the distance from a query point to the nearest point on a
    Great Circle Arc
    Args:
        lat1 (float): Latitude of start of arc in degrees
        lon1 (float): Longtiude of start of arc in degrees
        lat2 (float): Latitude of end of arc in degrees
        lon2 (float): Longtiude of end of arc in degrees
        plat (float): Latitude of query point in degrees
        plon (float): Longitude of query point in degrees
    Returns:
        float: The distance in kilometres from the query point to the great circle arc
    """
    tlat,tlon = NearestPointOnGC(alat,alon,blat,blon,plat,plon) #Nearest pt on GC
    abdist = Haversine(alat,alon,blat,blon)  #Length of arc
    atdist = Haversine(alat,alon,tlat,tlon)  #Distance arc start to nearest pt
    tbdist = Haversine(tlat,tlon,blat,blon)  #Distance arc end to nearest pt
    #If the nearest point T on the Great Circle lies within the arc, then the
    #length of the arc is approximately equal to the distance from T to each end
    #of the arc, accounting for floating-point errors
    PRECISION = 1e-3 #km 
    #We set the precision to a relatively high value because super-accuracy is not
    #to needed here and a failure to catch this can lead to vast under-estimates
    #of distance
    if np.abs(abdist-atdist-tbdist)<PRECISION: #Nearest point was on the arc
        return Haversine(tlat,tlon,plat,plon)
    #Okay, the nearest point wasn't on the arc, so the nearest point is one of the
    #ends points of the arc
    apdist = Haversine(alat,alon,plat,plon)
    bpdist = Haversine(blat,blon,plat,plon)
    return min(apdist,bpdist)


def Distance3dPointTo3dPolygon(lat,lon,geom):
    """
    Calculate the closest distance between a latitude-longitude query point
    and a `shapely` polygon defined by latitude-longitude points using only 
    spherical mathematics.

    Args:
        lat (float): Latitude of query point in degrees
        lon (float): Longitude of query point in degrees
        geom (`shapely` geometry): A `shapely` geometry whose points are in 
            latitude-longitude space
    Returns:
        dist (float): The minimum distance in kilometres between the polygon 
            and the query point.
    """
    if geom.type.values[0] == 'Polygon':
        dist = math.inf
#        xy = features[0]['geometry'][0].exterior.xy
        xy = mapping(geom)["features"][0]["geometry"]["coordinates"][0]
        #Polygons are closed rings, so the first-last pair is automagically delt with
        for p1, p2 in Pairwise(zip(*xy)):
            dist = min(dist,DistanceToGCArc(p1[1],p1[0],p2[1],p2[0],lat,lon))
    elif geom.type == 'MultiPolygon':
        dist = min(*[Distance3dPointTo3dPolygon(lat,lon,part) for part in geom])
    return dist

###############################################################################

# Jo wrote this

def Distance3dPointTo3dCoords(lat, lon, polylat, polylon):
    """
    Calculate the closest distance between a latitude-longitude query point
    and a `shapely` polygon defined by latitude-longitude points using only 
    spherical mathematics, modified from Distance3dPointTo3dPolygon.
    
    Args:
        lat (float): Latitude of query point in degrees
        lon (float): Longitude of query point in degrees
        polylat (array): Latitude borders polygon 
        polylon (array): Longitude borders polygon 
        geom (`shapely` geometry): A `shapely` geometry whose points are in 
            latitude-longitude space
    Returns:
        dist (float): The minimum distance in kilometres between the polygon 
            and the query point.
    """
    
    dist = math.inf
    for i in range(len(polylat)):
        if i > len(polylat) - 2:
            break
        dist = min(dist, DistanceToGCArc(polylat[i], polylon[i], 
                    polylat[i+1], polylon[i+1], lat, lon))
    
    return dist


def saa_distance_shapely(lat, lon, geom):
    """
    Wrapper to calculate distance to shapely geometry (defined in
    latitude-longitude space) for a series of points.
    
    Args:
        lat (float): Latitude of query point in degrees
        lon (float): Longitude of query point in degrees
        geom (`shapely` geometry): A `shapely` geometry whose points are in 
            latitude-longitude space
    Returns:
        dist (float): The minimum distance in kilometres between the polygon 
            and the query point.
    """
    dists = []
    for i in range(len(lat)):
        dist = Distance3dPointTo3dPolygon(lat[i], lon[i], geom)
        dists.append(dist)
    return dists


def saa_distance(lat, lon, poly_lat, poly_lon):
    """
    Wrapper to calculate distance to polygon (defined in
    latitude-longitude space) for a series of points.
    
    Args:
        lat (float): Latitude of query point in degrees
        lon (float): Longitude of query point in degrees
        polylat (array): Latitude borders polygon 
        polylon (array): Longitude borders polygon 
    Returns:
        dist (float): The minimum distance in kilometres between the polygon 
            and the query point.
    """
    dists = []
    for i in range(len(lat)):
        dist = Distance3dPointTo3dCoords(lat[i], lon[i], poly_lat, poly_lon)
        dists.append(dist)
    return dists
