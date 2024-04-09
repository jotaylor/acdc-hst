import os
import yaml
import getpass
from sqlalchemy.pool import NullPool
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

def load_connection(dbname, echo=False):
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

    cwd = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(cwd, "settings.yaml"), "r") as f:
        settings = yaml.load(f, Loader=yaml.SafeLoader)
        dbsettings = settings["dbsettings"][dbname]

    if dbsettings["dbtype"] == "sqlite":
        connection_string = f"sqlite:///{dbsettings['loc']}"
    elif dbsettings["dbtype"] == "mysql":
        pswd = getpass.getpass(f"Enter password for {dbsettings['dbname']}: ")
        connection_string = f"mysql+pymysql://{dbsettings['user']}:{pswd}@{dbsettings['host']}:{dbsettings['port']}/{dbsettings['dbname']}"
    engine = create_engine(connection_string, echo=echo, poolclass=NullPool)
    Session = sessionmaker(bind=engine)
    session = Session()

    return session, engine

