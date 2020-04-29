import yaml
import os
import glob
from astropy.io import fits
from sqlalchemy import Table

from connect_db import load_connection
from schema import Darks
from calculate_dark import measure_darkrate, dark_edges

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)

def populate_darks(files, connection_string=SETTINGS["connection_string"], 
                   db_table=Darks):
    session, base, engine = load_connection(connection_string)
    darks_table = Table(db_table.__tablename__, base.metadata, autoload=True)
    
    for i,item in enumerate(files):
        with fits.open(item) as hdulist:
            hdr0 = hdulist[0].header
            hdr1 = hdulist[1].header
        itemname = os.path.basename(item)
        dark_df = measure_darkrate(item)
        
        for j in range(len(dark_df["darkrate"])):
            dark_row = dark_df.iloc[j]
            darks_data = [
                {"filename": itemname,
                 "segment": hdr0["segment"],
                 "rootname": hdr0["rootname"],
                 "expstart": hdr1["expstart"],
                 "exptime": hdr1["exptime"],
                 "hva": hdr1["hvlevela"],
                 "hvb": hdr1["hvlevelb"],
                 "fileloc": item,
                 "region": dark_row["region"],
                 "xcorr_min": dark_row["xcorr_min"],
                 "xcorr_max": dark_row["xcorr_max"],
                 "ycorr_min": dark_row["ycorr_min"],
                 "ycorr_max": dark_row["ycorr_max"],
                 "latitude": dark_row["latitude"],
                 "longitude": dark_row["longitude"],
                 "saa_flag": dark_row["saa_flag"],
                 "darkrate": dark_row["darkrate"],
                 "time": dark_row["time"],
                 "unfiltered_pha_counts": dark_row["unfiltered_pha_counts"]}]
        darks_table.insert().execute(darks_data)
    print("Added {} rows to table Darks".format(len(files)))

if __name__ == "__main__":
    files = glob.glob(os.path.join(SETTINGS["data_dir"], "*corrtag*fits"))
    populate_darks(files)

