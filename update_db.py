from timeit import default_timer as timer
import yaml
import os
import glob
import numpy as np
from astropy.io import fits
from sqlalchemy import Table
from sqlalchemy.ext.declarative import declarative_base

from connect_db import load_connection
from schema import Solar, Darks
from calculate_dark import measure_darkrate, dark_edges, parse_solar_files

TESTING = True
TIMING = True

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)

def populate_darks(files, connection_string=SETTINGS["connection_string"], 
                   db_table=Darks):
    session, engine = load_connection(connection_string)
    base = declarative_base(engine)
    darks_table = Table(db_table.__tablename__, base.metadata, autoload=True)
    for i,item in enumerate(files):
        with fits.open(item) as hdulist:
            hdr0 = hdulist[0].header
            hdr1 = hdulist[1].header
        itemname = os.path.basename(item)

        dark_df = measure_darkrate(item)
        try:
            solar_q = session.query(Solar.time, Solar.flux)\
                        .filter(Solar.time >= min(dark_df["time"])-2)\
                        .filter(Solar.time <= max(dark_df["time"])+2)
            results = solar_q.all()
            solar_time = [x.time for x in results]
            solar_flux = [x.flux for x in results]
            solar_flux_interp = np.interp(dark_df["time"], solar_time, solar_flux)
        except:
            solar_flux_interp = np.full(len(dark_df["darkrate"]), -999)

        for j in range(len(dark_df["darkrate"])):
            start = timer()
            dark_row = dark_df.iloc[j]
            dark_data = [
                {"filename": itemname,
                 "segment": hdr0["segment"],
                 "rootname": hdr0["rootname"],
                 "expstart": hdr1["expstart"],
                 "exptime": hdr1["exptime"],
                 "hva": hdr1["hvlevela"],
                 "hvb": hdr1["hvlevelb"],
                 "fileloc": item,
                 "region": dark_row["region"],
                 "xcorr_min": int(dark_row["xcorr_min"]),
                 "xcorr_max": int(dark_row["xcorr_max"]),
                 "ycorr_min": int(dark_row["ycorr_min"]),
                 "ycorr_max": int(dark_row["ycorr_max"]),
                 "latitude": dark_row["latitude"],
                 "longitude": dark_row["longitude"],
                 "saa_flag": int(dark_row["saa_flag"]),
                 "darkrate": dark_row["darkrate"],
                 "time": dark_row["time"],
                 "unfiltered_pha_counts": dark_row["unfiltered_pha_counts"],
                 "solar_flux": solar_flux_interp[j]}]
            if TESTING is True:
                darks_table.insert().execute(dark_data)
                end = timer()
                if TIMING is True:
                    print("One insert took {} seconds".format(end-start))
                break
            darks_table.insert().execute(dark_data)
            end = timer()
            if TIMING is True:
                print("One insert took {} seconds".format(end-start))
    print("Updated table Darks")

def populate_solar(files, connection_string=SETTINGS["connection_string"],
                   db_table=Solar):
    session, engine = load_connection(connection_string)
    base = declarative_base(engine)
    solar_table = Table(db_table.__tablename__, base.metadata, autoload=True)

    for item in files:
        time, flux = parse_solar_files(item)
        for i in range(len(time)): 
            start = timer()
            solar_data = [
                {"time": time[i],
                 "flux": flux[i]}]
            if TESTING is True:
                solar_table.insert().execute(solar_data)
                end = timer()
                if TIMING is True:
                    print("One insert took {} seconds".format(end-start))
                break
            solar_table.insert().execute(solar_data)
            end = timer()
            if TIMING is True:
                print("One insert took {} seconds".format(end-start))
    print("Updated table Solar")
    
#def add_solar_to_darks():
        

if __name__ == "__main__":
    solar_files = glob.glob(os.path.join(SETTINGS["solar_dir"]))
    populate_solar(solar_files)
    files = glob.glob(os.path.join(SETTINGS["dark_dir"], "lbo2dkkuq_*corrtag*fits"))
    populate_darks(files)
