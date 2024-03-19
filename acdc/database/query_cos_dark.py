import yaml
from sqlalchemy.orm import load_only
import pandas as pd

from .connect_db import load_connection
from .schema import Solar, Darks

TESTING = False

def all_darks(dbname="cos_dark"):
    """
    Query the Darks table to get darks for all time for all PHAs.

    Args:
        dbname (str): The location of the SQLite database, with full path, e.g.
            /path/to/cos_dark.db 
            If in the current directory, do not include . or ./ 
    Returns:
        df (:obj:`pandas.DataFrame`): Query results.
    """
    
	# Connect to database
    session, engine = load_connection(dbname)

    # Define columns to return from database
    cols = ["expstart", "solar_flux", "latitude", "longitude", "segment", "hv", "region", "saa_distance"]
    cols += [f"dark_pha{x}" for x in range(0,32)]
    if TESTING is True: 
        cols = ["expstart", "latitude", "longitude", "solar_flux", "dark_pha11"]

    cols_attr = [getattr(Darks, f) for f in cols]
    # Execute SELECt query with WHERE statements
    query = session.query(Darks).options(load_only(*cols_attr))
#                .filter(Darks.region == "inner")\
#                .filter(Darks.segment == "FUVA")
    if TESTING is True:
        query = query.filter(Darks.filename == "le0a3mcuq_corrtag_a.fits")
    results = query.all()

    # Convert query results (list of attributes) to pandas dataframe
    d = {col: [getattr(x, col) for x in results] for col in cols}
    df = pd.DataFrame(d)
    
    return df

def files_by_mjd(mjdstart, mjdend, segment="FUVA", hv=167, morecols=[],
                 dbname="cos_dark"):
    
	# Connect to database
    session, engine = load_connection(dbname)

    cols = ["fileloc"]

    if hv == "*" or hv == None:
        query = session.query(Darks.fileloc.distinct())\
                .filter(Darks.expstart >= mjdstart)\
                .filter(Darks.expstart < mjdend)\
                .filter(Darks.segment == segment)\
                .order_by(Darks.expstart)
    else:
        query = session.query(Darks.fileloc.distinct())\
                .filter(Darks.expstart >= mjdstart)\
                .filter(Darks.expstart < mjdend)\
                .filter(Darks.segment == segment)\
                .filter(Darks.hv == hv)\
                .order_by(Darks.expstart)
    results = query.all()
    
    d = {cols[i]: [x[i] for x in results] for i in range(len(cols))}
    df = pd.DataFrame(d)

    return df

def counts_by_mjd(mjdstart, mjdend, morecols=[],
                  dbname="cos_dark"):
    
	# Connect to database
    session, engine = load_connection(dbname)
    cols = ["expstart", "solar_flux", "latitude", "longitude", "segment", "hv", "region", "saa_distance"]
    cols += [f"dark_pha{x}" for x in range(0,32)]
    cols += morecols
    cols_attr = [getattr(Darks, f) for f in cols]

    query = session.query(Darks).options(load_only(*cols_attr))\
                .filter(Darks.expstart >= mjdstart)\
                .filter(Darks.expstart < mjdend)
    results = query.all()

    d = {col: [getattr(x, col) for x in results] for col in cols}
    df = pd.DataFrame(d)

    return df

