import yaml
from sqlalchemy.orm import load_only
from sqlalchemy.sql import text
import pandas as pd

from connect_db import load_connection
import darkevents_schema

DBNAME = "dark_events"
COLUMNS = ["xcorr", "ycorr", "pha", "mjd", "hv", "segment", "filename", "proposid"]
TIMING = False

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


def all_rows(hvtable, returncols=COLUMNS):
    """
    Query the Darks table to get darks for all time for all PHAs.

    Args:
    Returns:
        df (:obj:`pandas.DataFrame`): Query results.
    """
    
    # Connect to database
    session, engine = load_connection(DBNAME)
    tablename = f"DarkEventsHv{hvtable}"
    events_table = getattr(darkevents_schema, tablename)

    # Execute query
    query = session.query(events_table).options(load_only(*returncols))
    results = query.all()

    df = sql_to_df(results, returncols)
    
    return df


def sql_query(hvtable, sqlquery):
    session, engine = load_connection(DBNAME)
    tablename = f"DarkEventsHv{hvtable}"
    events_table = getattr(darkevents_schema, tablename)

    query = session.query(events_table).from_statement(text(sqlquery))
    results = query.all()

    df = sql_to_df(results, COLUMNS)
    
    return df


def equals_query(hvtable, returncols=COLUMNS, **kwargs):
    session, engine = load_connection(DBNAME)
    tablename = f"DarkEventsHv{hvtable}"
    events_table = getattr(darkevents_schema, tablename)

    condition = {col: True for col in COLUMNS}
    if len(kwargs) != 0:
        for col in COLUMNS:
            if col in kwargs:
                colattr = getattr(events_table, col)
                condition[col] = colattr == kwargs[col]

    query = session.query(events_table).options(load_only(*returncols))\
                .filter(condition["xcorr"])\
                .filter(condition["ycorr"])\
                .filter(condition["pha"])\
                .filter(condition["mjd"])\
                .filter(condition["hv"])\
                .filter(condition["segment"])\
                .filter(condition["filename"])\
                .filter(condition["proposid"])
    results = query.all()
    
    df = sql_to_df(results, returncols)
   
    return df 


def range_query(hvtable, returncols=COLUMNS, segment="*", mjdstart=0, mjdend=99999, 
                x0=0, x1=16385, y0=0, y1=1025, pha0=0, pha1=32, hv="*", 
                proposid="*", filename="*"):
    session, engine = load_connection(DBNAME)
    tablename = f"DarkEventsHv{hvtable}"
    events_table = getattr(darkevents_schema, tablename)

    special = ["segment", "hv", "proposid", "filename"]
    condition = {}
    for col in special:
        val = locals()[col] # Get the value of the same-named variable
        if val == "*":
            condition[col] = True
        else:
            colattr = getattr(events_table, col)
            condition[col] = colattr == val

    query = session.query(events_table).options(load_only(*returncols))\
                .filter(condition["segment"])\
                .filter(condition["hv"])\
                .filter(events_table.mjd > mjdstart)\
                .filter(events_table.mjd < mjdend)\
                .filter(events_table.xcorr > x0)\
                .filter(events_table.xcorr < x1)\
                .filter(events_table.ycorr > y0)\
                .filter(events_table.ycorr < y1)\
                .filter(events_table.pha > pha0)\
                .filter(events_table.pha < pha1)\
                .filter(condition["filename"])\
                .filter(condition["proposid"])
    results = query.all()
    
    df = sql_to_df(results, returncols)
   
    return df 

