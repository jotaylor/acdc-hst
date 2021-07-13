import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
from scipy.stats import norm
from sqlalchemy.orm import load_only
import pandas as pd

from connect_db import load_connection
from darkevents_schema import DarkEvents

HV = {"FUVA": [163, 167, 169, 171, 173, 178],
      "FUVB": [163, 167, 169, 175]}
ALL_PHAS = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]

def query_data(segment, mjdstart, mjdend, x0, x1, y0, y1, hv=163, getall=False, dbname="cos_dark"):

    # Connect to database
    with open("settings.yaml", "r") as f:
        settings = yaml.load(f)
        dbsettings = settings["dbsettings"][dbname]
    session, engine = load_connection(dbsettings)
    
    # Define columns to return from database
    cols = ["pha", "hv"]
    
    # Execute query
    if getall is True:
        query = session.query(DarkEvents).options(load_only(*cols))\
                .filter(DarkEvents.hv == hv)
    else:
        query = session.query(DarkEvents).options(load_only(*cols))\
                .filter(DarkEvents.segment == segment)\
                .filter(DarkEvents.hv == hv)\
                .filter(DarkEvents.mjd > mjdstart)\
                .filter(DarkEvents.mjd < mjdend)\
                .filter(DarkEvents.xcorr > x0)\
                .filter(DarkEvents.xcorr < x1)\
                .filter(DarkEvents.ycorr > y0)\
                .filter(DarkEvents.ycorr < y1)
    results = query.all()

    # Convert query results (list of attributes) to pandas dataframe
    d = {col: [getattr(x, col) for x in results] for col in cols}
    df = pd.DataFrame(d)

    return df

def bin_data(df, segment):
    x = []
    y = []
    hvs_gain = defaultdict(list)
    for hv in HV[segment]:
        phas = df.loc[df["hv"] == hv].pha.values
        s_phas = list(set(phas))
        present = [x for x in ALL_PHAS if x in s_phas]
        if set(present) != set(ALL_PHAS):
            print(f"HV={hv} does not have counts at all PHAs")
            continue
        hist, bins = np.histogram(phas, bins=np.arange(1,32))
        mean,std = norm.fit(hist)
        modal_gain = mean
        hvs_gain[hv].append(modal_gain)
    for i in range(len(HV[segment])):
        if i == 0 or HV[segment][i] not in hvs_gain:
            continue
        previous_hvs = np.arange(i)
        for j in previous_hvs:
            hv_diff = HV[segment][i] - HV[segment][j]
            gain_diff = hvs_gain[HV[segment][i]] - hvs_gain[HV[segment][j]]
            relation = gain_diff / hv_diff
            starting_gain = hvs_gain[HV[segment][j]]
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
