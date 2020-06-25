import yaml
from sqlalchemy.orm import load_only
import pandas as pd

from connect_db import load_connection
from schema import Solar, Darks

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)

def query_darks(connection_string=SETTINGS["connection_string"]):
    # Connect to database
    session, engine = load_connection(connection_string)

    cols = ["expstart", "region"]
    cols += [f"dark_pha{x}" for x in range(0,32)]
#    cols = ["dark_pha10"]
    query = session.query(Darks).options(load_only(*cols))\
                .filter(Darks.region == "inner")
#                .filter(Darks.filename == "le0a3mcuq_corrtag_a.fits")
#                .filter(Darks.id == 9)
    results = query.all()

    d = {col: [getattr(x, col) for x in results] for col in cols}
    df = pd.DataFrame(d)
    
    return df

