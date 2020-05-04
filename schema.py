from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column
from sqlalchemy import Float, ForeignKey, Integer, String

Base = declarative_base()

class Darks(Base):
    """Object-relational mapping for the Darks table."""

    __tablename__ = "darks"
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
    ycorr_min = Column(Integer())
    ycorr_max = Column(Integer())
    latitude = Column(Float())
    longitude = Column(Float())
    saa_flag = Column(Integer())
    darkrate = Column(Float())
    time = Column(Float())
    unfiltered_pha_counts = Column(Float())
    solar_flux = Column(Float())

class Solar(Base):
    """Object-relational mapping for the Solar table."""

    __tablename__ = "solar"
    id = Column(Integer, primary_key=True, nullable=False)
    time = Column(Float())
    flux = Column(Float())

