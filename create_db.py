import yaml
from connect_db import load_connection
from schema import Base, Darks, Solar

with open("settings.yaml", "r") as f:
    SETTINGS = yaml.load(f)

def create_db(connection_string=SETTINGS["connection_string"]):
    session, engine = load_connection(SETTINGS["connection_string"])
    # It's important that the Base from schema.py be used (from the import)
    Base.metadata.create_all(engine)
    print("Created database cos_dark.db")

if __name__ == "__main__":
    create_db()
