from timeit import default_timer as timer
import yaml
import os
import glob
from astropy.io import fits
from sqlalchemy import Table
from sqlalchemy.ext.declarative import declarative_base

from connect_db import load_connection
from schema import Solar, Darks
from calculate_dark import measure_darkrate, dark_edges, parse_solar_files

TESTING = True

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

        start = timer()
        dark_df = measure_darkrate(item)
        end = timer()
        print("measure_darkrate() took {} seconds".format(end-start))
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
                 "unfiltered_pha_counts": dark_row["unfiltered_pha_counts"]}]
            if TESTING is True:
                darks_table.insert().execute(dark_data)
                end = timer()
                print("One insert took {} seconds".format(end-start))
                break
            darks_table.insert().execute(dark_data)
            end = timer()
            print("One insert took {} seconds".format(end-start))
    print("Updated table Darks")

def populate_solar(files, connection_string=SETTINGS["connection_string"],
                   db_table=Solar):
    session, engine = load_connection(connection_string)
    base = declarative_base(engine)
    solar_table = Table(db_table.__tablename__, base.metadata, autoload=True)

    for item in files:
        start = timer()
        date, flux = parse_solar_files(item)
        end = timer()
        print("parse_solar_files() took {} seconds".format(end-start))
        for i in range(len(date)): 
            start = timer()
            solar_data = [
                {"time": date[i],
                 "flux": flux[i]}]
            if TESTING is True:
                solar_table.insert().execute(solar_data)
                end = timer()
                print("One insert took {} seconds".format(end-start))
                break
            solar_table.insert().execute(solar_data)
            end = timer()
            print("One insert took {} seconds".format(end-start))
    print("Updated table Solar")
            

if __name__ == "__main__":
    files = glob.glob(os.path.join(SETTINGS["dark_dir"], "lbo2dkkuq_*corrtag*fits"))
    populate_darks(files)
    solar_files = glob.glob(os.path.join(SETTINGS["solar_dir"]))
    populate_solar(solar_files)
