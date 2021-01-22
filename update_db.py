import datetime
import yaml
import os
import glob
import numpy as np
from timeit import default_timer as timer
from astropy.io import fits
from sqlalchemy import Table
from sqlalchemy.ext.declarative import declarative_base

import create_db
from connect_db import load_connection
from schema import Solar, Darks
from calculate_dark import measure_darkrate, parse_solar_files, get_1291_box
from within_saa import get_saa_poly
from saa_distance import Distance3dPointTo3dCoords

# For testing purposes only. If TESTING = True, only one value is recorded per
# input dataset to save time. If TIMING = True, recorded runtime for each insert
# is written to STDOUT.
TESTING = False
TIMING = False

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)
    DBNAME = SETTINGS["dbname"]

def populate_solar(files, dbname=DBNAME, db_table=Solar):
    """
    Populate the Solar database table. This can be used to add new files
    or ingest all data from scratch.
    In order to properly get interpolated solar flux values for dark
    observations, Solar table *must* be populated before the Darks table.
    
    Args:
        files (array-like): Names of solar flux files to parse.
        dbname (str): The location of the SQLite database, with full path, e.g.
            /path/to/cos_dark.db 
            If in the current directory, do not include . or ./ 
        db_table (:obj:`sqlalchemy.Table`): Name of table to update, 
            defaults to Solar.
    """
    
    # Connect to database.
    session, engine = load_connection(dbname)
    base = declarative_base(engine)
    solar_table = Table(db_table.__tablename__, base.metadata, autoload=True)

    for fileno,item in enumerate(files):
        loop0 = timer()    
        all_solar_rows = []

        time, flux = parse_solar_files(item)
        for i in range(len(time)): 
            start = timer()
            solar_data = {"time": time[i],
                          "flux": flux[i]}
            all_solar_rows.append(solar_data)

        insert0 = timer()
        solar_table.insert().execute(all_solar_rows)
        insert1 = timer()
        print("File {}/{}".format(fileno+1, len(files))) 

        if TIMING is True:
            print("One insert took {:.3g} seconds".format(insert1-insert0))
            print("One file loop took {:.3g} seconds".format(insert1-loop0))
        if TESTING is True:
            print("Updated table Solar")
            return

    print("Updated table Solar")

def populate_darks(files, dbname=DBNAME, db_table=Darks):
    """
    Populate the Darks database table. This can be used to add new files
    or ingest all data from scratch.
    In order to properly get interpolated solar flux values for dark
    observations, Solar table *must* be populated before the Darks table.
    
    Args:
        files (array-like): Names of COS dark files to ingest.
        dbname (str): The location of the SQLite database, with full path, e.g.
            /path/to/cos_dark.db 
            If in the current directory, do not include . or ./ 
        db_table (:obj:`sqlalchemy.Table`): Name of table to update,
            defaults to Darks.
    """

    psa_1291 = get_1291_box()

    # Connect to database.
    session, engine = load_connection(dbname)
    base = declarative_base(engine)
    darks_table = Table(db_table.__tablename__, base.metadata, autoload=True)
    
    # Handle distance to SAA
    saa_dists = {}
    saa_lat, saa_lon = get_saa_poly()

    for i,item in enumerate(files):
        loop0 = timer()
        all_dark_rows = []

        with fits.open(item) as hdulist:
            hdr0 = hdulist[0].header
            hdr1 = hdulist[1].header
        itemname = os.path.basename(item)
        info, dark = measure_darkrate(item, psa_1291)
        
        # Query the Solar table to get solar fluxes near the dark observation
        # dates. Interpolate to get appropriate values. If Solar table is not
        # properly up to date, insert -999 values.
        try:
            solar_q = session.query(Solar.time, Solar.flux)\
                        .filter(Solar.time >= min(info["time"])-2)\
                        .filter(Solar.time <= max(info["time"])+2)
            results = solar_q.all()
            solar_time = [x.time for x in results]
            solar_flux = [x.flux for x in results]
            solar_flux_interp = np.interp(info["time"], solar_time, solar_flux)
        except:
            solar_flux_interp = np.full(len(info["time"]), -999)

        # For each time sampling of the dark observation file, insert a series of
        # values into the Darks table row.
        # Loop over time indices (by default, 25s intervals)
        dark_data0 = {"filename": itemname,
           "segment": info["segment"],
           "exptime": info["exptime"],
           "hv": hdr1["HVLEVEL{}".format(info["segment"][-1])],
           "fileloc": item}
        for j in range(len(info["time"])):
            dark_data1 = dark_data0.copy()
            lat = round(info["latitude"][j], 6)
            lon = round(info["longitude"][j], 6)
            if (lat, lon) not in saa_dists:
                dist = Distance3dPointTo3dCoords(lat, lon, saa_lat, saa_lon)
                saa_dists[(lat, lon)] = dist
            else:
                dist = saa_dists[(lat, lon)]
            dark_data1["latitude"] = lat
            dark_data1["longitude"] = lon
            dark_data1["saa_distance"] = dist
            dark_data1["solar_flux"] = solar_flux_interp[j]
            dark_data1["expstart"] = round(info["time"][j], 6)
            # Loop over regions.
            for region in dark.keys():
                dark_data = dark_data1.copy()
                dark_data["region"] = region
                dark_data["region_area"] = int(dark[region]["region_area"])
                dark_data["xcorr_min"] = int(dark[region]["xcorr_min"])
                dark_data["xcorr_max"] = int(dark[region]["xcorr_max"])
                dark_data["ycorr_min"] = int(dark[region]["ycorr_min"])
                dark_data["ycorr_max"] = int(dark[region]["ycorr_max"])
                dark_data["region_area"] = dark[region]["region_area"]
                for k in range(len(info["pha"])):
                    pha_num = "dark_pha{}".format(info["pha"][k])
                    dark_data[pha_num] = int(dark[region]["darks"][j][k])
                all_dark_rows.append(dark_data)
        insert0 = timer()
        darks_table.insert().execute(all_dark_rows)
        insert1 = timer()
        print("File {}/{}".format(i+1, len(files))) 
        
        if TIMING is True:
            print("One insert took {:.3g} seconds".format(insert1-insert0))
            print("One file loop took {:.3g} seconds".format(insert1-loop0))
        if TESTING is True:
            print("Updated table Darks")
            return

    print("Updated table Darks")

if __name__ == "__main__":
    if not os.path.exists("cos_dark.db"):
        create_db.create_db()
    start = datetime.datetime.now()
    print("Start time: {}".format(start))
    solar_files = glob.glob(os.path.join(SETTINGS["solar_dir"]))
    populate_solar(solar_files)
    files = glob.glob(os.path.join(SETTINGS["dark_dir"]))
    populate_darks(files)
    end = datetime.datetime.now()
    print("End time: {}".format(end))
    print("Total time: {}".format(end-start))
