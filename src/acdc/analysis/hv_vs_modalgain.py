import os
import yaml
import numpy as np
import matplotlib.pyplot as plt
dirname = os.path.dirname(__file__)
stylesheet = os.path.join(dirname, "niceplot.mplstyle")
plt.style.use(stylesheet)
from scipy.stats import norm
from sqlalchemy.orm import load_only
from timeit import default_timer as timer
import pandas as pd
from collections import OrderedDict,defaultdict

from acdc.database.connect_db import load_connection
from acdc.database import darkevents_schema
from acdc.utils.utils import sql_to_df

HV = {"FUVA": [163, 167, 169, 171, 173, 178],
      "FUVB": [163, 167, 169, 175]}
ALL_PHAS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]
DBNAME = "dark_events"
TIMING = True

def query_data(segment="*", mjdstart=0, mjdend=99999, x0=0, x1=16385, 
               y0=0, y1=1025, hv=163, getall=True):

    # Connect to database
    session, engine = load_connection("dark_events")
    
    # Define columns to return from database
    cols = ["pha", "hv", "segment"]

    # Get Tables for all HVs
#    possible_hvs = [163, 167, 169, 171, 173, 175, 178]
    possible_hvs = [171, 173]
    data = {}
    for hv in possible_hvs:
        tablename = f"DarkEventsHv{hv}"
        events_table = getattr(darkevents_schema, tablename)
        # Execute query
        query0 = timer()
        if getall is True:
            query = session.query(events_table).options(load_only(*cols))
        else:
            query = session.query(events_table).options(load_only(*cols))\
                    .filter(events_table.segment == segment)\
                    .filter(events_table.hv == hv)\
                    .filter(events_table.mjd > mjdstart)\
                    .filter(events_table.mjd < mjdend)\
                    .filter(events_table.xcorr > x0)\
                    .filter(events_table.xcorr < x1)\
                    .filter(events_table.ycorr > y0)\
                    .filter(events_table.ycorr < y1)
        results = query.all()
        query1 = timer()

        # Convert query results (list of attributes) to pandas dataframe
        df = sql_to_df(results, cols)
        data[hv] = df
        processing1 = timer()
        if TIMING is True:
            print(f"Query for HV={hv} took {query1-query0:.3g} seconds")
            print(f"Data formatting for HV={hv} took {processing1-query1:.3g} seconds")

    return data

def bin_data(data, segment, interactive=False):
    x = []
    y = []
    hvs_gain = OrderedDict()
    for hv in HV[segment]:
        df = data[hv] 
        phas = df.pha.values
        s_phas = list(set(phas))
        present = [x for x in ALL_PHAS if x in s_phas]
        if set(present) != set(ALL_PHAS):
            print(f"HV={hv} does not have counts at all PHAs")
            continue
        if interactive is True:
            plt.hist(phas, bins=np.arange(0,32))
            plt.show()
            x0 = int(input("Enter starting X for normal fit: "))
            x1 = int(input("Enter ending X for normal fit: "))
            plt.clf()
        else:
            x0 = 5
            x1 = 15
        subset = np.where((phas > x0) & (phas < x1))
        mean,std = norm.fit(phas[subset])
        modal_gain = mean
        hvs_gain[hv] = modal_gain
    hvs_keys = list(hvs_gain.keys())
    for i,hv in enumerate(hvs_keys):
        if i == 0:
            continue
        previous_hvs = np.arange(i)
        for j in previous_hvs:
            hv_diff = hv - hvs_keys[j] 
            gain_diff = hvs_gain[hv] - hvs_gain[hvs_keys[j]]
            relation = gain_diff / hv_diff
            starting_gain = hvs_gain[hvs_keys[j]]
            x.append(starting_gain)
            y.append(relation)
    
    return x,y

def plot_data(x, y, segment):
    fig, ax = plt.subplots(figsize=(24,12))
    ax.plot(x, y, "k,")
    ax.set_title("f{segment}")
    ax.set_xlabel("Starting Modal Gain")
    ax.set_ylabel("Pulse Height Bins per HV Step")
    out = f"{segment}_pha_hv_relation.png"
    plt.savefig(out, bbox_inches="tight")
    print(f"Wrote {out}")
