import os
from sqlalchemy.pool import NullPool
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

def load_connection(dbname, echo=False):
    """
    Return session, base, and engine objects for connecting to the database.

    Args:
        dbname (str): The location of the SQLite database, with full path, e.g.
            /path/to/cos_dark.db 
            If in the current directory, do not include . or ./ 
            This will get plugged into the full connection string of form:
            `dialect+driver://username:password@host:port/database`
        echo (Bool): If True, log all statements to STDOUT. Defaults to False.

    Returns:
        session (:obj:`Session()` instance): Provides a holding zone 
            for all objects loaded or associated with the database.
        engine (:obj:`engine` object): Provides a source of database 
            connectivity and behavior.
    """

    connection_string = os.path.join("sqlite:///", dbname)
    engine = create_engine(connection_string, echo=echo, poolclass=NullPool)
    Session = sessionmaker(bind=engine)
    session = Session()

    return session, engine

