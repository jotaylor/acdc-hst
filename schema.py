from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column
from sqlalchemy import Float, ForeignKey, Integer, String

Base = declarative_base()

class Darks(Base):
    """Object-relational mapping for the Darks table."""

    __tablename__ = "darks"
    id = Column(Integer, primary_key=True, nullable=False)
    filename = Column(String(30))
    rootname = Column(String(10))
    segment = Column(String(5))
    expstart = Column(Float())
    exptime = Column(Float())
    hva = Column(Integer())
    hvb = Column(Integer())
    latitude = Column(Float())
    longitude = Column(Float())
    darkrate = Column(Float())
    time = Column(Float())
    region = Column(String(10))
    xcorr_min = Column(Integer())
    xcorr_max = Column(Integer())
    ycorr_min = Column(Integer())
    ycorr_max = Column(Integer())
    saa_flag = Column(Integer())
    unfiltered_pha_counts = Column(Float())
    solar_flux = Column(Float())
    fileloc = Column(String(200))

class Solar(Base):
    """Object-relational mapping for the Solar table."""

    __tablename__ = "solar"
    id = Column(Integer, primary_key=True, nullable=False)
    time = Column(Float())
    flux = Column(Float())

