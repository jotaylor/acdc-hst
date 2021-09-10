from sqlalchemy import Index
import yaml

from connect_db import load_connection
from darkevents_schema import DarkEvents

# Connect to database
with open("settings.yaml", "r") as f:
    settings = yaml.load(f, Loader=yaml.SafeLoader)
    dbsettings = settings["dbsettings"]["dark_events"]
session, engine = load_connection(dbsettings)

hv_index = Index("hv_index", DarkEvents.hv)

hv_index.create(bind=engine)
