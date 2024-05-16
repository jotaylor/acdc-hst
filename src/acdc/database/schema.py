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
    expstart = Column(Float())
    exptime = Column(Float())
    hv = Column(Integer())
    latitude = Column(Float())
    longitude = Column(Float())
    saa_distance= Column(Float())
    region = Column(String(15))
    region_area = Column(Integer())
    xcorr_min = Column(Integer())
    xcorr_max = Column(Integer())
    ycorr_min = Column(Integer())
    ycorr_max = Column(Integer())
    solar_flux = Column(Float())
    fileloc = Column(String(200))
    dark_pha0 = Column(Integer())
    dark_pha1 = Column(Integer())
    dark_pha2 = Column(Integer())
    dark_pha3 = Column(Integer())
    dark_pha4 = Column(Integer())
    dark_pha5 = Column(Integer())
    dark_pha6 = Column(Integer())
    dark_pha7 = Column(Integer())
    dark_pha8 = Column(Integer())
    dark_pha9 = Column(Integer())
    dark_pha10 = Column(Integer())
    dark_pha11 = Column(Integer())
    dark_pha12 = Column(Integer())
    dark_pha13 = Column(Integer())
    dark_pha14 = Column(Integer())
    dark_pha15 = Column(Integer())
    dark_pha16 = Column(Integer())
    dark_pha17 = Column(Integer())
    dark_pha18 = Column(Integer())
    dark_pha19 = Column(Integer())
    dark_pha20 = Column(Integer())
    dark_pha21 = Column(Integer())
    dark_pha22 = Column(Integer())
    dark_pha23 = Column(Integer())
    dark_pha24 = Column(Integer())
    dark_pha25 = Column(Integer())
    dark_pha26 = Column(Integer())
    dark_pha27 = Column(Integer())
    dark_pha28 = Column(Integer())
    dark_pha29 = Column(Integer())
    dark_pha30 = Column(Integer())
    dark_pha31 = Column(Integer())

class Solar(Base):
    """Object-relational mapping for the Solar table."""

    __tablename__ = "solar"
    id = Column(Integer, primary_key=True, nullable=False)
    time = Column(Float())
    flux = Column(Float())

