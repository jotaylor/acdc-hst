import argparse
import yaml
from connect_db import load_connection

DBNAME = "cos_dark.db" 

def create_sqlite_db(dbname=DBNAME):
    """
    Create the database, with name specified from `settings.yaml`.

    Args:
        connection_string (str): Connection string of the form:
            `dialect+driver://username:password@host:port/database`
    """

    if dbname == "cos_dark.db":
        from schema import Base, Darks, Solar
    with open("settings.yaml", "r") as f:
        settings = yaml.load(f, Loader=yaml.SafeLoader)
        dbsettings = settings["dbsettings"][dbname]
    session, engine = load_connection(dbsettings)
    # It's important that the Base from schema.py be used (from the import)
    Base.metadata.create_all(engine)
    print(f"Created database {dbname}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--db", default=DBNAME,
                        help="Name of database to create")
    args = parser.parse_args()
    create_sqlite_db(args.db)
