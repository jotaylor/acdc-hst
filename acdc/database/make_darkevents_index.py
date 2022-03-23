from sqlalchemy import Index
import yaml

from .connect_db import load_connection
from .darkevents_schema import DarkEvents

# Connect to database
session, engine = load_connection("dark_events")

hv_index = Index("hv_index", DarkEvents.hv)

hv_index.create(bind=engine)
