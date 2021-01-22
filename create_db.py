import yaml
from connect_db import load_connection
from schema import Base, Darks, Solar

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)
    DBNAME = SETTINGS["dbname"]

def create_db(dbname=DBNAME):
    """
    Create the database, with name specified from `settings.yaml`.

    Args:
        connection_string (str): Connection string of the form:
            `dialect+driver://username:password@host:port/database`
    """

    session, engine = load_connection(dbname)
    # It's important that the Base from schema.py be used (from the import)
    Base.metadata.create_all(engine)
    print("Created database cos_dark.db")

if __name__ == "__main__":
    create_db()
