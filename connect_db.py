import os
import getpass
from sqlalchemy.pool import NullPool
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

def load_connection(dbsettings, echo=False):
    """
    Return session, base, and engine objects for connecting to the database.

    Args:
        dbsettings (dict): Dictionary listing database type and all required
            connection parameters.
            Params will get plugged into the full connection string of form:
            `dialect+driver://username:password@host:port/database`
        echo (Bool): If True, log all statements to STDOUT. Defaults to False.

    Returns:
        session (:obj:`Session()` instance): Provides a holding zone 
            for all objects loaded or associated with the database.
        engine (:obj:`engine` object): Provides a source of database 
            connectivity and behavior.
    """

    if dbsettings["dbtype"] == "sqlite":
        connection_string = os.path.join("sqlite://", dbname)
    elif dbsettings["dbtype"] == "mysql":
        pswd = getpass.getpass(f"Enter password for {dbsettings['dbname']}: ")
        connection_string = f"mysql+pymysql://{dbsettings['user']}:{pswd}@{dbsettings['host']}:{dbsettings['port']}/{dbsettings['dbname']}"
    engine = create_engine(connection_string, echo=echo, poolclass=NullPool)
    Session = sessionmaker(bind=engine)
    session = Session()

    return session, engine

