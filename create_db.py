import yaml
from connect_db import load_connection
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column
from sqlalchemy import Float, ForeignKey, Integer, String

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)
session, base, engine = load_connection(SETTINGS["connection_string"])

class Darks(base):
    """ORM for the header table"""

    __tablename__ = 'darks'
    id = Column(Integer, primary_key=True, nullable=False)
    filename = Column(String(30))
    segment = Column(String(5))
    rootname = Column(String(10))
    expstart = Column(Float())
    exptime = Column(Float())
    hva = Column(Integer())
    hvb = Column(Integer())
    fileloc = Column(String(200))
    region = Column(String(10))
    xcorr_min = Column(Integer())
    xcorr_max = Column(Integer())
    ycorr_max = Column(Integer())
    ycorr_max = Column(Integer())
    latitude = Column(Float())
    longitude = Column(Float())
    saa_flag = Column(Integer())
    darkrate = Column(Float())
    time = Column(Float())
    unfiltered_pha_counts = Column(Float())

def create_db():
	base.metadata.create_all()
	print("Created SQLite database ./cos_dark.db")

if __name__ == "__main__":
    create_db()
