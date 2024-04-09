from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy import Column
from sqlalchemy import Float, ForeignKey, Integer, String
from sqlalchemy.types import DECIMAL as Decimal

Base = declarative_base()

class DarkEvents(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkevents"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())

class DarkEventsHv163(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshv163"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())

class DarkEventsHv167(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshv167"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())

class DarkEventsHv169(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshv169"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())

class DarkEventsHv171(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshv171"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())

class DarkEventsHv173(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshv173"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())

class DarkEventsHv175(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshv175"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())

class DarkEventsHv178(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshv178"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())


class DarkEventsHvOther(Base):
    """Object-relational mapping for the DarkEvents table."""

    __tablename__ = "darkeventshvother"
    id = Column(Integer, primary_key=True, nullable=False)
    xcorr = Column(Decimal(12,5))
    ycorr = Column(Decimal(12,5))
    pha = Column(Integer())
    mjd = Column(Decimal(precision='20,11'))
    hv = Column(Integer())
    segment = Column(String(5))
    filename = Column(String(30))
    proposid = Column(Integer())
