import argparse
import yaml

from acdc.database.connect_db import load_connection

def create_db(dbname="cos_dark"):
    """
    Create the database, with name specified from `settings.yaml`.

    Args:
        dbname (str): Name of database to create with empty tables.
            Currently cos_dark and dark_events are supported. 
    """

    if dbname == "cos_dark":
        from schema import Base, Darks, Solar
    elif dbname == "dark_events":
        from darkevents_schema import Base, DarkEvents
    session, engine = load_connection(dbname)
    # It's important that the Base from appropriate be schema be used, imported above
    Base.metadata.create_all(engine)
    print(f"Created database & tables for {dbname}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--db", default="cos_dark",
                        help="Name of database to create")
    args = parser.parse_args()
    create_db(args.db)
