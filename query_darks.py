import yaml
from sqlalchemy.orm import load_only
import pandas as pd

from connect_db import load_connection
from schema import Solar, Darks

TESTING = False

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)

def all_darks(connection_string=SETTINGS["connection_string"]):
    """
    Query the Darks table to get darks for all time for all PHAs.

    Args:
        connection_string (str): Connection string to DB of the form:
            `dialect+driver://username:password@host:port/database`
    Returns:
        df (:obj:`pandas.DataFrame`): Query results.
    """
    
	# Connect to database
    session, engine = load_connection(connection_string)

    # Define columns to return from database
    cols = ["expstart", "solar_flux", "latitude", "longitude"]
    cols += [f"dark_pha{x}" for x in range(0,32)]
    if TESTING is True: 
        cols = ["expstart", "latitude", "longitude", "solar_flux", "dark_pha11"]

    # Execute SELECt query with WHERE statements
    query = session.query(Darks).options(load_only(*cols))\
                .filter(Darks.region == "inner")\
                .filter(Darks.segment == "FUVA")
    if TESTING is True:
        query = query.filter(Darks.filename == "le0a3mcuq_corrtag_a.fits")
    results = query.all()

    # Convert query results (list of attributes) to pandas dataframe
    d = {col: [getattr(x, col) for x in results] for col in cols}
    df = pd.DataFrame(d)
    
    return df

def files_by_mjd(mjdstart, mjdend, morecols=[],
                 connection_string=SETTINGS["connection_string"]):
    
	# Connect to database
    session, engine = load_connection(connection_string)

    cols = ["fileloc", "segment"]

    query = session.query(Darks.fileloc.distinct(), Darks.segment)\
                .filter(Darks.expstart >= mjdstart)\
                .filter(Darks.expstart < mjdend)
    results = query.all()

    d = {cols[i]: [x[i] for x in results] for i in range(len(cols))}
    df = pd.DataFrame(d)

    return df

def counts_by_mjd(mjdstart, mjdend, morecols=[],
                  connection_string=SETTINGS["connection_string"]):
    
	# Connect to database
    session, engine = load_connection(connection_string)

    cols = ["segment", "exptime", "hv"]
    cols += [f"dark_pha{x}" for x in range(0,32)]
    cols += morecols

    query = session.query(Darks).options(load_only(*cols))\
                .filter(Darks.expstart >= mjdstart)\
                .filter(Darks.expstart < mjdend)
    results = query.all()

    d = {col: [getattr(x, col) for x in results] for col in cols}
    df = pd.DataFrame(d)

    return df

