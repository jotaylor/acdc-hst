import argparse
import yaml
from connect_db import load_connection

def create_sqlite_db(dbname="cos_dark"):
    """
    Create the database, with name specified from `settings.yaml`.

    Args:
        connection_string (str): Connection string of the form:
            `dialect+driver://username:password@host:port/database`
    """

    if dbname == "cos_dark":
        from schema import Base, Darks, Solar
    elif dbname == "dark_events":
        from darkevents_schema import Base, DarkEvents
    with open("settings.yaml", "r") as f:
        settings = yaml.load(f, Loader=yaml.SafeLoader)
        dbsettings = settings["dbsettings"][dbname]
    session, engine = load_connection(dbsettings)
    # It's important that the Base from appropriate be schema be used, imported above
    Base.metadata.create_all(engine)
    print(f"Created database/tables for {dbname}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--db", default="cos_dark",
                        help="Name of database to create")
    args = parser.parse_args()
    create_sqlite_db(args.db)
