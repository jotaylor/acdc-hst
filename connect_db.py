from sqlalchemy.pool import NullPool
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

def load_connection(connection_string, echo=False):
    """
    Return session, base, and engine objects for connecting to the database.

    Args:
        connection_string (str): The connection string to connect to the 
            database. The connection string should take the form:
            `dialect+driver://username:password@host:port/database`
        echo (Bool): If True, log all statements to STDOUT. Defaults to False.

    Returns:
        session (:obj:`Session()` instance): Provides a holding zone 
            for all objects loaded or associated with the database.
        engine (:obj:`engine` object): Provides a source of database 
            connectivity and behavior.
    """

    engine = create_engine(connection_string, echo=echo, poolclass=NullPool)
    Session = sessionmaker(bind=engine)
    session = Session()

    return session, engine

