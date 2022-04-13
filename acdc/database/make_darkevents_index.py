from sqlalchemy import Index

from acdc.database.connect_db import load_connection
from acdc.database.darkevents_schema import DarkEvents

# Connect to database
session, engine = load_connection("dark_events")

hv_index = Index("hv_index", DarkEvents.hv)

hv_index.create(bind=engine)
